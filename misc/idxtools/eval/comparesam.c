#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <getopt.h>
#include <unistd.h>
#include <pthread.h>
#include <signal.h>
#include <errno.h> 
#include <sys/resource.h>
#include <sys/time.h>


#include <string>
#include <iostream>
using namespace std;

//maximum number of mappings for a particular read
#define NUM_MAPPINGS 300  

//SAM flags
#define F_UNMAPPED 0x04
#define F_SECONDARY 0x100
#define F_SUPPLEMENTARY 0x800

/*
by default, only exact matches are taken to be correct (i.e. chr and the exact position)
define OVERLAP_BASED_EVAL to allow an overlap
*/
//#define OVERLAP_BASED_EVAL 1

/* 
by default, only primary mappings are considered. define CONSIDER_SUPPLEMENTARY
to consider supplementary reads set CONSIDER_SUPPLEMENTARY
*/
//#define CONSIDER_SUPPLEMENTARY 1

/*
by default, secondary mappings are ignored. define CONSIDER_SECONDARY to consider supplementary.
however, this option works only if CONSIDER_SUPPLEMENTARY is also set
*/
//#define CONSIDER_SECONDARY 1 


//used for debugging
#define DEBUG_LOOP 1

//SAM entry
typedef struct{
    string rid; //read name
    string tid; //reference sequence name
    int32_t target_start; //mapping coordinate start
    int32_t target_end; //mapping coordinate end
    int32_t mapq; //MAPQ
    int32_t score; //AS
    int32_t flag; //sam flag
    char strand; //+ or - strand
}alignment_t;

//sam file
typedef struct{

    //sam file statistics
    int32_t num_entries=0; //number of sam entries
	int32_t mapped_reads=0; //total number of reads that map (primary)
	int32_t mappings=0; //total number of alignments 
    int32_t unmapped=0;
    int32_t secondary=0;
    int32_t supplementary=0;

    //temporary data structure thats hold all the mappings for the current read
    alignment_t mapping_of_reads[NUM_MAPPINGS]; 
    int32_t mapping_of_reads_n=0;

} sam_stat_t;


//global statistic of the sam comparison - between samfile 'a' and 'b'
typedef struct{

    //number of reads that fall under each scenario
    int32_t unmapped_in_both=0; //the read is unmapped in both 'a' and 'b'
    int32_t same=0;             //correct - the read maps to the same location (locations overlaps if overlap based evaluation is set)
    int32_t mismatches=0;       //in correct
    int32_t only_in_a=0;        //read is only mapped in  'a'
    int32_t only_in_b=0;        //read is only mapped in  'b'
    int32_t pri_a_to_supp_b=0;  //primary mapping in 'a' is a supplementary mapping in 'b'
    int32_t pri_a_to_sec_b=0;   //primary mapping in 'b' is a secondary mappings in 'b'

    //file pointers for tsv
    FILE* f_only_in_a;
    FILE* f_only_in_b;
    FILE* f_mismatches;

    //file pointers for bed
    FILE* f_only_in_a_bed;
    FILE* f_only_in_b_bed;
    FILE* f_mismatches_bed_a;
    FILE* f_mismatches_bed_b;

} compare_stat_t;






#define MALLOC_CHK(ret) malloc_chk((void*)ret, __func__, __FILE__, __LINE__ - 1)
#define F_CHK(ret, filename)                                                   \
    f_chk((void*)ret, __func__, __FILE__, __LINE__ - 1, filename);
#define NULL_CHK(ret) null_chk((void*)ret, __func__, __FILE__, __LINE__ - 1)
#define NEG_CHK(ret) neg_chk(ret, __func__, __FILE__, __LINE__ - 1)

static inline void malloc_chk(void* ret, const char* func, const char* file,
                              int line) {
    if (ret != NULL)
        return;
    fprintf(
        stderr,
        "[%s::ERROR]\033[1;31m Failed to allocate memory : "
        "%s.\033[0m\n[%s::DEBUG]\033[1;35m Error occured at %s:%d.\033[0m\n\n",
        func, strerror(errno), func, file, line);
    exit(EXIT_FAILURE);
}

static inline void f_chk(void* ret, const char* func, const char* file,
                         int line, const char* fopen_f) {
    if (ret != NULL)
        return;
    fprintf(
        stderr,
        "[%s::ERROR]\033[1;31m Failed to open %s : "
        "%s.\033[0m\n[%s::DEBUG]\033[1;35m Error occured at %s:%d.\033[0m\n\n",
        func, fopen_f, strerror(errno), func, file, line);
    exit(EXIT_FAILURE);
}

// Die on error. Print the error and exit if the return value of the previous function NULL
static inline void null_chk(void* ret, const char* func, const char* file,
                            int line) {
    if (ret != NULL)
        return;
    fprintf(stderr,
            "[%s::ERROR]\033[1;31m %s.\033[0m\n[%s::DEBUG]\033[1;35m Error "
            "occured at %s:%d.\033[0m\n\n",
            func, strerror(errno), func, file, line);
    exit(EXIT_FAILURE);
}

// Die on error. Print the error and exit if the return value of the previous function is -1
static inline void neg_chk(int ret, const char* func, const char* file,
                           int line) {
    if (ret >= 0)
        return;
    fprintf(stderr,
            "[%s::ERROR]\033[1;31m %s.\033[0m\n[%s::DEBUG]\033[1;35m Error "
            "occured at %s:%d.\033[0m\n\n",
            func, strerror(errno), func, file, line);
    exit(EXIT_FAILURE);
}




void insert_alignment(sam_stat_t *read, alignment_t *sam_entry){

#ifndef CONSIDER_SUPPLEMENTARY
    if((sam_entry->flag & (F_UNMAPPED|F_SECONDARY|F_SUPPLEMENTARY))==0){ //only primary
#else
    #ifndef CONSIDER_SECONDARY
        if((sam_entry->flag & (F_UNMAPPED|F_SECONDARY))==0){ //primary+supplementary
    #else
        if((sam_entry->flag & (F_UNMAPPED))==0){ //primary+supplementary+secondary
    #endif
#endif
        if(read->mapping_of_reads_n>=NUM_MAPPINGS){
            fprintf(stderr,"Too many mapping for read %s\n",sam_entry->rid.c_str());
            fprintf(stderr,"Flag %d\n",sam_entry->flag);
            read->mapping_of_reads_n--;
        }           
        read->mapping_of_reads[read->mapping_of_reads_n]=*sam_entry;
        read->mapping_of_reads_n++;

    }   

}


//get the end coordinate of a mapping
//from https://www.biostars.org/p/17891/
static inline int32_t cigar_to_aln(char *cigar,int32_t pos){
    char *c = cigar;
    int32_t pos_end=pos;
    while(*c!=0){
        char* p2;
        if(!isdigit(*c)){
            fprintf(stderr,"bad cigar string1: %s\n",cigar);
            exit(EXIT_FAILURE);
        }
        int32_t len=strtol(c,&p2,10);
        if(len<=0){
            fprintf(stderr,"bad cigar string2: %s\n",cigar);
            exit(EXIT_FAILURE);
        }
        switch(*p2){
            case 'M':
            {
                pos_end += len; 
                break;
            }
            case 'D':
            {
                pos_end += len;
                break;
            }
            // case 'I':
            // {
            //    n_gap += len; 
            // }

            default:
            {
                break;
            }
        }
        c=p2+1;
    }
    return pos_end;
}


void parse_sam_entry(alignment_t *sam_entry, char *buffer){
        //read name
        char *pch = strtok (buffer,"\t\r\n"); assert(pch!=NULL);
        sam_entry->rid = pch;

        //flag
        pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
        sam_entry->flag = atoi(pch);

        //if mapped
        if(!(sam_entry->flag & 0x04)){        
            sam_entry->strand = (sam_entry->flag&16)? '-' : '+';

            //RNAME
            pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
            sam_entry->tid = pch;
            
            //POS
            pch = strtok (NULL,"\t\r\n"); assert(pch!=NULL);
            sam_entry->target_start = atoi(pch);
            
            //MAPQ
            pch = strtok (NULL,"\t\r\n");  assert(pch!=NULL);
            sam_entry->mapq = atoi(pch);
            
            //CIGAR
            pch = strtok (NULL,"\t\r\n");  assert(pch!=NULL);
            sam_entry->target_end = cigar_to_aln(pch,sam_entry->target_start);

            //RNEXT
            pch = strtok (NULL,"\t\r\n");  assert(pch!=NULL);
            
            //PNEXT
            pch = strtok (NULL,"\t\r\n");  assert(pch!=NULL);
            
            //TLEN
            pch = strtok (NULL,"\t\r\n");  assert(pch!=NULL);
            
            //SEQ
            pch = strtok (NULL,"\t\r\n");  assert(pch!=NULL);
            
            //QUAL
            pch = strtok (NULL,"\t\r\n");  assert(pch!=NULL);
        
            //NM
            pch = strtok (NULL,"\t\r\n");  assert(pch!=NULL);
            
            //ms
            pch = strtok (NULL,"\t\r\n");  assert(pch!=NULL);
            
            //alignment score
            pch = strtok (NULL,"\t\r\n");  assert(pch!=NULL);
            int32_t as;
            int32_t ret=sscanf(pch,"AS:i:%d",&as);
            assert(ret>0);
            sam_entry->score=as;
        }
}

//mapping in 'a' and 'b' are exactly the same?
static inline int8_t is_correct_exact(alignment_t a, alignment_t b){
    if(a.tid==b.tid && a.target_start==b.target_start && a.strand==b.strand){
        return 1;
    }
    else{
        return 0;
    }

}

#ifdef OVERLAP_BASED_EVAL 
//adapted from paftools
/*Suppose the reported mapping coordinate in samfile 'a' overlap with the 
mapping coordinate in samfile 'b'  as follows:

mapping in 'a':   --------------------
mapping in 'b':           ----------------------
                  |<- l1 ->|<-- o -->|<-- l2 -->|
Let r=o/(l1+o+l2). The reported mapping is considered correct if r>0.1 by default.
*/
#define ovlp_ratio 0.1
int8_t is_correct_overlap(alignment_t a, alignment_t b)
{
    if (a.tid != b.tid || a.strand != b.strand) return 0;
    float o, l;
    if (a.target_start < b.target_start) {
        if (a.target_end <= b.target_start) return 0;
        o = (a.target_end < b.target_end? a.target_end : b.target_end) - b.target_start;
        l = (a.target_end > b.target_end? a.target_end : b.target_end) - a.target_start;
    } else {
        if (b.target_end <= a.target_start) return 0;
        o = (a.target_end < b.target_end? a.target_end : b.target_end) - a.target_start;
        l = (a.target_end > b.target_end? a.target_end : b.target_end) - b.target_start;
    }
    return o/l > ovlp_ratio? 1 : 0;
}
#endif

//compare 
void compare_alnread(compare_stat_t *compare, sam_stat_t *reada,sam_stat_t *readb){
    if(reada->mapping_of_reads_n==0 && readb->mapping_of_reads_n==0){
        compare->unmapped_in_both++;
        return;
    }
    //only in b
    if(readb->mapping_of_reads_n>0 && reada->mapping_of_reads_n==0){
        compare->only_in_b++;    
        alignment_t b = readb->mapping_of_reads[0];

        //tsv - 1 based index
        fprintf(compare->f_only_in_b,"%s\t%s\t%d\t%d\t%d\n",b.rid.c_str(),
        b.tid.c_str(),b.target_start,b.mapq, b.score);

        //bed - 0 based index for star coordinate, 1 based index for end coordinate
        fprintf(compare->f_only_in_b_bed,"%s\t%d\t%d\t%s\t%d\t%c\t%d\n",b.tid.c_str(),b.target_start-1,
        b.target_end,b.rid.c_str(),b.score,b.strand,b.mapq);

        return;
    }
    //only in a
    if(reada->mapping_of_reads_n>0 && readb->mapping_of_reads_n==0){
        compare->only_in_a++;
        alignment_t a = reada->mapping_of_reads[0];

        //tsv - 1 based index    
        fprintf(compare->f_only_in_a,"%s\t%s\t%d\t%d\t%d\n",a.rid.c_str(),
        a.tid.c_str(),a.target_start,a.mapq, a.score);

        //bed - 0 based index for star coordinate, 1 based index for end coordinate
        fprintf(compare->f_only_in_a_bed,"%s\t%d\t%d\t%s\t%d\t%c\t%d\n",a.tid.c_str(),a.target_start-1,
        a.target_end,a.rid.c_str(),a.score,a.strand,a.mapq);

        return;
    }

    //compare the mappings for the read in 'a' with the mappings for the same read in 'b'
    int32_t flag=0;
    for(int32_t i=0;i < reada->mapping_of_reads_n ;i++){
        for(int32_t j=0; j < readb->mapping_of_reads_n ;j++){
                alignment_t a = reada->mapping_of_reads[i];
                alignment_t b = readb->mapping_of_reads[j];
                if(a.rid!=b.rid){ //read names should match
                    cout << a.rid << '\t' << b.rid <<endl;
                }
                assert(a.rid==b.rid);
        //primary only mode        
        #ifndef CONSIDER_SUPPLEMENTARY        
            #ifndef OVERLAP_BASED_EVAL 
                //exact mapping or not?
                if(is_correct_exact(a,b)){
                    flag++; 
                }
            #else
                //coordinates overlap or not?
                if(is_correct_overlap(a,b)){
                    flag++; 
                }
            #endif
        #else
                //in 'a' we always take the primary only
                if((a.flag&(F_SUPPLEMENTARY|F_SECONDARY))==0){
    
                    #ifndef OVERLAP_BASED_EVAL 
                        //same mapping
                        if(is_correct_exact(a,b)){
                            flag++;
                            if(b.flag&F_SECONDARY){
                                compare->pri_a_to_sec_b++;
                            }
                            else if(b.flag&F_SUPPLEMENTARY){
                                compare->pri_a_to_supp_b++;
                            }
                            break;
                        }
                    #else
                        if(is_correct_overlap(a,b)){
                            flag++; 
                            if(b.flag&F_SECONDARY){
                                compare->pri_a_to_sec_b++;
                            }
                            else if(b.flag&F_SUPPLEMENTARY){
                                compare->pri_a_to_supp_b++;
                            }
                            break;
                        }
                    #endif


                }
        #endif    
        }
    }  

    assert(flag==0 || flag ==1);
    if(flag){
        compare->same++; //correct
    }
    else{
        compare->mismatches++; //incorrect
        alignment_t a = reada->mapping_of_reads[0];
        alignment_t b = readb->mapping_of_reads[0];

        //1 based index for tsv
        fprintf(compare->f_mismatches,"%s\t%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\n",a.rid.c_str(),
        a.tid.c_str(),a.target_start,a.mapq, a.score,b.tid.c_str(),b.target_start,b.mapq, b.score);

        //bed - 0 index start ,1 based end
        fprintf(compare->f_mismatches_bed_a,"%s\t%d\t%d\t%s\t%d\t%c\t%d\n",a.tid.c_str(),a.target_start-1,
        a.target_end,a.rid.c_str(),a.score,a.strand,a.mapq);
        fprintf(compare->f_mismatches_bed_b,"%s\t%d\t%d\t%s\t%d\t%c\t%d\n",b.tid.c_str(),b.target_start-1,
        b.target_end,b.rid.c_str(),b.score,b.strand,b.mapq);

    }

}

void print_sam_stat(char * filename, sam_stat_t* a){
    printf("File %s\n"
    "Number of entries\t%d\n"
    "Total number of alignments\t%d\n"
    "Number of mapped reads\t%d\n"
    "Unmapped\t%d\n"
    "Secondary\t%d\n"
    "Supplementary\t%d\n"
    ,filename,a->num_entries,a->mappings,a->mapped_reads,a->unmapped, a->secondary, a->supplementary);
}

void print_compare_stat(compare_stat_t *compare){
    printf("\nComparison between a and b\n"
    "Unmapped in both\t%d\n"
    "Same\t%d\n"
    "Mismatches\t%d\n"
    "Only in A\t%d\n"
    "Only in B\t%d\n"
    ,compare->unmapped_in_both,compare->same, compare->mismatches, compare->only_in_a, compare->only_in_b);
#ifdef CONSIDER_SUPPLEMENTARY
    printf("Primary in 'a' is a supplementary in 'b'\t%d\n",compare->pri_a_to_supp_b);
#endif
#ifdef CONSIDER_SECONDARY
    printf("Primary in 'a' is a secondary in 'b'\t%d\n",compare->pri_a_to_sec_b);
#endif

}

//global sam stat update
void update_sam_stat(sam_stat_t* a,alignment_t curr_a){
    a->num_entries++;
    if(!(curr_a.flag & F_UNMAPPED)){
        a->mappings++;
        if(curr_a.flag & F_SECONDARY){
            a->secondary++;
        }
        else if(curr_a.flag & F_SUPPLEMENTARY){
            a->supplementary++;
        }
        else{
            a->mapped_reads++;
        }
    }
    else{
        a->unmapped++;
    }
            
}

/*
Let the 'a' and 'b' be two samfiles for the same set of reads, 
but mapped with different mappers (or same mapper with different options)

compare sam will compare 'a' and 'b' and give mappings statistics such as the number of reads  :
- unmapped in both
- correct : the read maps to the same location (locations overlaps if overlap based evaluation is set)
- incorrect 
- read is only mapped in  'a'
- read is only mapped in  'b'
- primary mapping in 'a' is a supplementary mapping in 'b'
- primary mapping in 'b' is a secondary mappings in 'b'

compare sam will also output (as tsv files and bed file) reads that :
-mismatch between 'a' and 'b' 
-unique to 'a'
-unique to 'b'

functionality : 

samfile 'a' and 'b' are sequentially read while loading all the mappings for a 
particular read name at a time. 
IMPORTANT : the samfiles 'a' and 'b' should have the reads in the same order and
the multiple mappings for a given read ID should be adjacently located.

For each loaded read, we compare the mappings between 'a' and 'b'
In 'a' always the primary mapping is considered despite the value of 
CONSIDER_SUPPLEMENTARY and CONSIDER_SECONDARY flags.
In 'b' supplementary and secondary mappings can also be considered 
if the flags CONSIDER_SUPPLEMENTARY and CONSIDER_SECONDARY are set.

The comparison statistics are updated during the comparison and 
the entries will be written by the tsv and bed file if required.

Finally we print the comparison statistics.

*/


int main(int argc, char* argv[]){
	
    optind=1;

    if (argc-optind < 2) {
        fprintf(
            stderr,
            "Usage: %s [OPTIONS] a.sam b.sam\n",
            argv[0]);
        exit(EXIT_FAILURE);
    }

    FILE* samfile_a= fopen(argv[optind],"r");
    F_CHK(samfile_a,argv[optind]);
    FILE* samfile_b= fopen(argv[optind+1],"r");
    F_CHK(samfile_b,argv[optind+1]);


    /****************************sam A**********************************/
    //buffers for getline
    size_t bufferSize_a = 4096;
    char *buffer_a = (char *)malloc(sizeof(char)*(bufferSize_a)); 
    MALLOC_CHK(buffer_a);
    int readlinebytes_a=1;

#ifdef DEBUG_LOOP
    int32_t num_entries_a=0;
	int32_t mapped_reads_a=0;
	int32_t mappings_a=0;
#endif 

    alignment_t prev_a;
    prev_a.rid = "";
    alignment_t curr_a;
    curr_a.rid = "";

    sam_stat_t a;

    /****************************sam b********************************/
    //buffers for getline
    size_t bufferSize_b = 4096;
    char *buffer_b = (char *)malloc(sizeof(char)*(bufferSize_b)); 
    MALLOC_CHK(buffer_b);
    int readlinebytes_b=1;
    //char *pch_b=NULL;

#ifdef DEBUG_LOOP
    int32_t num_entries_b=0;
    int32_t mapped_reads_b=0;
    int32_t mappings_b=0;
#endif

    alignment_t prev_b;
    prev_b.rid = "";
    alignment_t curr_b;
    curr_b.rid = "";

    sam_stat_t b;

    compare_stat_t compare;
    compare.f_only_in_a=fopen("only_in_a.tsv","w");
    fprintf(compare.f_only_in_a,"readID\ttarget\ttarget_pos\tmapq\tAS\n");
    compare.f_only_in_b=fopen("only_in_b.tsv","w");
    fprintf(compare.f_only_in_b,"readID\ttarget\ttarget_pos\tmapq\tAS\n");
    compare.f_mismatches=fopen("mismatches.tsv","w");
    fprintf(compare.f_mismatches,"readID\ttarget_a\ttarget_pos_a\tmapq_a\tAS_a\ttarget_b\ttarget_pos_b\tmapq_b\tAS_b\n");

    compare.f_only_in_a_bed=fopen("only_in_a.bed","w");
    compare.f_only_in_b_bed=fopen("only_in_b.bed","w");
    compare.f_mismatches_bed_a=fopen("mismatches_a.bed","w");
    compare.f_mismatches_bed_b=fopen("mismatches_b.bed","w");



    int8_t state=0;

    /** Read all mappings for the cuurrent read in SAM A*/
	while(1){
        readlinebytes_a=getline(&buffer_a, &bufferSize_a, samfile_a); 
        if(readlinebytes_a == -1){
            assert(getline(&buffer_b, &bufferSize_b, samfile_b)==-1);
            if(state==1){
                compare_alnread(&compare,&a,&b);  
            }
            break;
        } 

        //ignore header lines
        if(buffer_a[0]=='@' || buffer_a[0]=='\n' || buffer_a[0]=='\r'){
            continue;
        }
        
        parse_sam_entry(&curr_a, buffer_a);
        update_sam_stat(&a,curr_a);

    #ifdef DEBUG_LOOP
        if(!(curr_a.flag & 0x04)){
            mappings_a++;
        }
        num_entries_a++;
    #endif

        
        if(curr_a.rid==prev_a.rid ){
            //cout <<  "enter a " << curr_a.rid<< endl;
            insert_alignment(&a, &curr_a);  
        }

        else{
            
        #ifdef DEBUG_LOOP
            if(!(curr_a.flag & 0x04)){
                mapped_reads_a++;
            }
        #endif

            int8_t break_flag = 0;

            /** Now read the entries for the same read in SAM B */   
            while(1){
                readlinebytes_b=getline(&buffer_b, &bufferSize_b, samfile_b); 
                if(readlinebytes_b == -1){
                    break_flag=-1;
                    break;
                } 

                //ignore header lines
                if(buffer_b[0]=='@' || buffer_b[0]=='\n' || buffer_b[0]=='\r'){
                    continue;
                }
                
                parse_sam_entry(&curr_b, buffer_b);
                update_sam_stat(&b,curr_b);
            #ifdef DEBUG_LOOP    
                if(!(curr_b.flag & 0x04)){
                    mappings_b++;
                }    
                num_entries_b++;
            #endif

                if(curr_b.rid==prev_b.rid ){
                    //cout <<  "enter b " << curr_b.rid<< endl;
                    insert_alignment(&b, &curr_b);
                    
                }
                else{
                    prev_b = curr_b;
                    break_flag=1;
                    break;
                }
                
                prev_b = curr_b;

            }

            //do the comparison for the mappings in the read    
            if(state==1){
                compare_alnread(&compare,&a,&b);  
                state=0;
            }
        

            a.mapping_of_reads_n=0;
            //cout <<  "enter a2 " << curr_a.rid<< endl;
            insert_alignment(&a,&curr_a);
            state=1;
           

            assert(break_flag!=0);
            if(break_flag==1){
            #ifdef DEBUG_LOOP      
                 if(!(curr_b.flag & 0x04)){
                     mapped_reads_b++;
                 }
            #endif
                b.mapping_of_reads_n=0;
                //cout <<  "enter  b2 " << curr_b.rid << endl;

                insert_alignment(&b,&curr_b);
            }

        }
        prev_a = curr_a;
    }

#ifdef DEBUG_LOOP
    fprintf(stderr,"[%s] Num entries in a %d, mappings_a %d, mapped reads in a %d, Num_entries in b %d, mappings_b %d, mapped reads in b %d\n",
            __func__, num_entries_a, mappings_a,mapped_reads_a,num_entries_b, mappings_b,mapped_reads_b);

    assert(a.num_entries==num_entries_a);
    assert(b.num_entries==num_entries_b);
    assert(a.mappings==mappings_a);
    assert(b.mappings==mappings_b);
    assert(a.mapped_reads==mapped_reads_a);
    assert(b.mapped_reads==mapped_reads_b);
#endif

    assert(a.num_entries == a.unmapped + a.mapped_reads + a.secondary + a.supplementary);
    assert(a.mappings==a.mapped_reads + a.secondary + a.supplementary);
    assert(b.num_entries == b.unmapped + b.mapped_reads + b.secondary + b.supplementary);
    assert(b.mappings==b.mapped_reads + b.secondary + b.supplementary);

    assert(compare.unmapped_in_both+compare.only_in_a==b.unmapped);
    assert(compare.unmapped_in_both+compare.only_in_b==a.unmapped);
    assert(compare.same+compare.mismatches+compare.only_in_a==a.mapped_reads);
    assert(compare.same+compare.mismatches+compare.only_in_b==b.mapped_reads);
    

    print_sam_stat(argv[optind],&a);
    print_sam_stat(argv[optind+1],&b);
    print_compare_stat(&compare);

    fclose(compare.f_mismatches_bed_a);
    fclose(compare.f_mismatches_bed_b);
    fclose(compare.f_only_in_a_bed);
    fclose(compare.f_only_in_b_bed);
    fclose(compare.f_mismatches);
    fclose(compare.f_only_in_a);
    fclose(compare.f_only_in_b);
    fclose(samfile_a);
    fclose(samfile_b);
    free(buffer_a);
    free(buffer_b);
	
	return 0;
}
