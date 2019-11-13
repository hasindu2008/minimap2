#include <assert.h>
#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_CHR 1024 //maximum number of chromosomes
#define MAX_CHR_NAME                                                           \
    256 //maximum number of characters in a chromosome name (with the comments)
#define OUTPUT_FILE_FORMAT "part%d.fa" //output fasta file name format
#define FASTA_STAT "stat_fasta.csv"    //csv to write the fasta stats
#define PARTITION_STAT                                                         \
    "stat_part%d.csv" //csv file name format to store stats for each partition

int8_t verbose = 1;          //verbosity level 0,1,2
int8_t print_fasta_stat = 0; //whether to output FASTA_STAT or not
int32_t NUM_PARTS = 0;       //number of partitions in the index

/*Die on error. Print the error and exit if the return value of the previous function NULL*/
#define errorCheckNULL(ret)                                                    \
    ({                                                                         \
        if (ret == NULL) {                                                     \
            fprintf(stderr, "Error at File %s line number %d : %s.\n",         \
                    __FILE__, __LINE__, strerror(errno));                      \
            exit(EXIT_FAILURE);                                                \
        }                                                                      \
    })

/*Die on error. Print the error and exit if the return value of the previous function is -1*/
#define errorCheck(ret)                                                        \
    ({                                                                         \
        if (ret < 0) {                                                         \
            fprintf(stderr, "Error at File %s line number %d : %s.\n",         \
                    __FILE__, __LINE__, strerror(errno));                      \
            exit(EXIT_FAILURE);                                                \
        }                                                                      \
    })

typedef struct {
    char chr_name[MAX_CHR_NAME]; //name of chromosome
    int64_t chr_len;             //length of the chromosome in bases
    int64_t chr_len_no_N; //length of the chromosome in bases (no N bases)
} chr_data_t;

//compare function for descending sort
int32_t cmpfunc(const void* a, const void* b) {
    chr_data_t* A = (chr_data_t*)a;
    chr_data_t* B = (chr_data_t*)b;
    return ((B->chr_len_no_N) - (A->chr_len_no_N));
}

//find the index of the minimum value in an array
int32_t minimum_index(int64_t* list, int32_t listsize) {
    int64_t minimum = INT64_MAX;
    int32_t min_index = -1;
    int32_t j;
    for (j = 0; j < listsize; j++) {
        if (list[j] < minimum) {
            min_index = j;
            minimum = list[j];
        }
    }
    assert(min_index >= 0);
    return min_index;
}

//find to which part a certain chromosome(chr_name) is belonging to
//in theory performance can be improved using a hash_table
//however a genome does not contain millions of chromosomes, so for the moment O(n^2) method
int32_t belong_to_which_part(char* chr_name, chr_data_t chr_per_part[][MAX_CHR],
                             int32_t* numchr_per_part) {
    int32_t i, j;
    for (i = 0; i < NUM_PARTS; i++) {
        for (j = 0; j < numchr_per_part[i]; j++) {
            if (strcmp(chr_per_part[i][j].chr_name, chr_name) == 0) {
                return i;
            }
        }
    }
    assert(0);
    return -1;
}

int main(int argc, char** argv) {
    if (argc != 3) {
        fprintf(stderr, "Usage : %s <reference.fa> <num_parts>\n", argv[0]);
        fprintf(stderr, "reference.fa - path to the fasta file containing the "
                        "reference genome");
        fprintf(stderr, "num_parts - number of partitions in the index");
        fprintf(stderr, "Example : $0 hg19.fa 4");
        exit(EXIT_FAILURE);
    }

    //open the fasta
    FILE* fasta = fopen(argv[1], "r");
    errorCheckNULL(fasta);

    //number of parts
    NUM_PARTS = atoi(argv[2]);
    if (NUM_PARTS < 2) {
        fprintf(stderr, "ERROR : Number of partitions should be equal or "
                        "greater than 2.\n");
        exit(EXIT_FAILURE);
    }

    //for getline
    size_t bufferSize = 1000;
    char* buffer = (char*)malloc(sizeof(char) * bufferSize);
    errorCheckNULL(buffer);
    size_t readlinebytes = 0;

    int32_t num_chr = 0; //keep track of number of chromosomes in the fasta
    chr_data_t theGenome[MAX_CHR]; //information for the whole genome

    //variables for current chromosome being processed
    int64_t chr_len = 0;
    int64_t chr_len_no_N = 0;
    char chr_name[MAX_CHR_NAME];

    /*********************** Read through the fast a and collect stats ******************************************/
    //can be made faster by reading an fasta.fai index if needed

    if (verbose >= 1) {
        fprintf(stderr, "INFO : Collecting chromosome stats.\n");
    }
    while (1) {
        readlinebytes = getline(&buffer, &bufferSize, fasta);
        if (readlinebytes == -1) { //EOF
            if (num_chr !=
                0) { //should save the data for the previous chromosome
                theGenome[num_chr - 1].chr_len = chr_len;
                theGenome[num_chr - 1].chr_len_no_N = chr_len_no_N;
                strcpy(theGenome[num_chr - 1].chr_name, chr_name);
                if (verbose >= 2) {
                    fprintf(stderr, "Parsed %s.\n",
                            theGenome[num_chr - 1].chr_name);
                }
            }
            break;
        }
        if (readlinebytes == 0) {
            fprintf(stderr, "ERROR : We read nothing. Something is wrong in "
                            "the fasta file.\n");
            exit(EXIT_FAILURE);
        }

        //a new chromosome
        if (buffer[0] == '>') {
            if (num_chr !=
                0) { //if not the first chromosome in the file, save the data of the previous chromosome that was processed
                theGenome[num_chr - 1].chr_len = chr_len;
                theGenome[num_chr - 1].chr_len_no_N = chr_len_no_N;
                strcpy(theGenome[num_chr - 1].chr_name, chr_name);
                if (verbose >= 2) {
                    fprintf(stderr, "Parsed %s.\n",
                            theGenome[num_chr - 1].chr_name);
                }
            }

            num_chr++;

            //reset lengths for the new chromosome
            chr_len = 0;
            chr_len_no_N = 0;

            if (readlinebytes - 1 > MAX_CHR_NAME) {
                fprintf(stderr,
                        "ERROR : Chromosome name (plus comments) is too larger "
                        "than the hard coded value, Increase "
                        "MAX_CHR_NAME in the code.\n");
                exit(EXIT_FAILURE);
            }
            strcpy(
                chr_name,
                &buffer
                    [1]); //copy the chromosome name except the ">". In theory only the chr name is needed and the other comments can be removed (for efficiency)
            if (chr_name[strlen(chr_name) - 1] == '\n' ||
                chr_name[strlen(chr_name) - 1] ==
                    '\r') { //unix and mac os newline style
                chr_name[strlen(chr_name) - 1] = '\0';
            } else {
                fprintf(stderr, "ERROR : New line character should be either "
                                "'\n' or '\r'.\n");
                exit(EXIT_FAILURE);
                ;
            }
            if (chr_name[strlen(chr_name) - 2] == '\r') { //windows new lines
                chr_name[strlen(chr_name) - 2] = '\0';
            }

            if (num_chr > MAX_CHR) {
                fprintf(stderr, "ERROR : Too many chromosomes than the hard "
                                "coded value, Increase MAX_CHR in the code.\n");
                exit(EXIT_FAILURE);
            }
        }

        //going through the chromosome
        else {
            int32_t i = 0;
            //go through all bases
            for (i = 0; i < readlinebytes; i++) {
                if (buffer[i] == 'A' || buffer[i] == 'C' || buffer[i] == 'G' ||
                    buffer[i] == 'T' || buffer[i] == 'a' || buffer[i] == 'c' ||
                    buffer[i] == 'g' || buffer[i] == 't') {
                    chr_len++;
                    chr_len_no_N++;
                } else if (buffer[i] == 'R' || buffer[i] == 'Y' ||
                           buffer[i] == 'K' || buffer[i] == 'M' ||
                           buffer[i] == 'S' || buffer[i] == 'W' ||
                           buffer[i] == 'B' || buffer[i] == 'D' ||
                           buffer[i] == 'H' || buffer[i] == 'V' ||
                           buffer[i] == 'r' || buffer[i] == 'y' ||
                           buffer[i] == 'k' || buffer[i] == 'm' ||
                           buffer[i] == 's' || buffer[i] == 'w' ||
                           buffer[i] == 'b' || buffer[i] == 'd' ||
                           buffer[i] == 'h' || buffer[i] == 'v') {
                    chr_len++;
                    chr_len_no_N++;
                } else if (buffer[i] == 'N' || buffer[i] == 'n') {
                    chr_len++;
                } else if (buffer[i] == '\n' || buffer[i] == '\r' ||
                           buffer[i] == '\0') {
                } else {
                    fprintf(stderr,
                            "WARNING : Invalid character found in %s : '%c'.\n",
                            chr_name, buffer[i]);
                    chr_len++;
                }
            }
        }
    }

    fclose(fasta);

    int64_t sum = 0;
    if (verbose >= 1) {
        fprintf(stderr, "INFO : %d chromosomes parsed.\n", num_chr);
    }

    if (print_fasta_stat) {
        //print stats
        FILE* stat = fopen(FASTA_STAT, "w");
        errorCheckNULL(stat);
        fprintf(stat, "chromosome_name,chromosome_length_with_N,chromosome_"
                      "length_without_N)\n");
        int32_t i;
        for (i = 0; i < num_chr; i++) {
            fprintf(stat, "\"%s\",%ld,%ld\n", theGenome[i].chr_name,
                    theGenome[i].chr_len, theGenome[i].chr_len_no_N);
            sum += theGenome[i].chr_len_no_N;
        }
        fclose(stat);
    }

    //sort based on lengths
    //long per_part = sum/NUM_PARTS;
    qsort(theGenome, num_chr, sizeof(chr_data_t), cmpfunc);

    //for(i=0;i<num_chr;i++){
    //    fprintf (stderr,"\"%s\",%ld,%ld\n",theGenome[i].chr_name, theGenome[i].chr_len, theGenome[i].chr_len_no_N);

    //}

    /*********************** Now do the partitioning******************************************/

    if (verbose >= 1) {
        fprintf(stderr, "INFO : Determining partitions.\n");
    }
    chr_data_t chr_per_part
        [NUM_PARTS]
        [MAX_CHR]; //per each part we need to store the chromosomes to be processed
    int32_t numchr_per_part[NUM_PARTS]; //number of chromosomes for each part
    int64_t length_per_part
        [NUM_PARTS]; //the total length of chromosomes (no N) for each part
    int32_t i;

    for (i = 0; i < NUM_PARTS; i++) {
        numchr_per_part[i] = 0;
        length_per_part[i] = 0;
    }

    int32_t j;
    int32_t min_index = 0;

    //go through all the chromosomes
    //add the current chromosome to the list with lowest sum
    for (i = 0; i < num_chr; i++) {
        min_index = minimum_index(length_per_part, NUM_PARTS);
        chr_per_part[min_index][numchr_per_part[min_index]] = theGenome[i];
        numchr_per_part[min_index]++;
        length_per_part[min_index] += theGenome[i].chr_len_no_N;
    }

    char filename[1024];

    //print the stats for each part
    for (i = 0; i < NUM_PARTS; i++) {
        sprintf(filename, PARTITION_STAT, i);
        FILE* stat = fopen(filename, "w");
        errorCheckNULL(stat);
        fprintf(stat, "chromosome_name,chromosome_length_with_N,chromosome_"
                      "length_without_N\n");
        for (j = 0; j < numchr_per_part[i]; j++) {
            fprintf(stat, "%s,%ld,%ld\n", chr_per_part[i][j].chr_name,
                    chr_per_part[i][j].chr_len,
                    chr_per_part[i][j].chr_len_no_N);
        }
        fclose(stat);
        if (verbose >= 1) {
            fprintf(stderr,
                    "INFO : Partition %d - %d chromosomes (%ld non-ambiguous "
                    "bases). See %s.\n",
                    i, numchr_per_part[i], length_per_part[i], filename);
        }
    }

    /*********************** Now write the fastas for each part******************************************/

    if (verbose >= 1) {
        fprintf(stderr, "INFO : Writing partitions.\n");
    }

    fasta = fopen(argv[1], "r");
    errorCheckNULL(fasta);

    FILE* outputs[NUM_PARTS];

    for (i = 0; i < NUM_PARTS; i++) {
        sprintf(filename, OUTPUT_FILE_FORMAT, i);
        outputs[i] = fopen(filename, "w");
        errorCheckNULL(outputs[i]);
    }

    int32_t part_index = -1;

    //read each line of the input file
    while (1) {
        readlinebytes = getline(&buffer, &bufferSize, fasta);
        if (readlinebytes == -1) {
            break;
        }
        if (readlinebytes == 0) {
            fprintf(stderr, "ERROR : We read nothing. Something is wrong in "
                            "the fasta file.\n");
            exit(EXIT_FAILURE);
        }

        if (buffer[0] == '>') {
            if (readlinebytes - 1 > MAX_CHR_NAME) {
                fprintf(stderr, "ERROR : Chromosome name too large, Increase "
                                "MAX_CHR_NAME.\n");
                exit(EXIT_FAILURE);
            }
            strcpy(chr_name, &buffer[1]);
            if (chr_name[strlen(chr_name) - 1] == '\n' ||
                chr_name[strlen(chr_name) - 1] ==
                    '\r') { //unix and mac newline style
                chr_name[strlen(chr_name) - 1] = '\0';
            } else {
                fprintf(stderr, "ERROR : New line character should be either "
                                "'\n' or '\r'.\n");
                exit(EXIT_FAILURE);
                ;
            }
            if (chr_name[strlen(chr_name) - 2] == '\r') { //windows new lines
                chr_name[strlen(chr_name) - 2] = '\0';
            }

            //find to which part this chromosome must be allocated to
            part_index =
                belong_to_which_part(chr_name, chr_per_part, numchr_per_part);
            fprintf(outputs[part_index], "%s", buffer);
            if (verbose >= 2) {
                fprintf(stderr, "Writing chromosome %s to partition %d.\n",
                        chr_name, part_index);
            }

        } else {
            //write to the correct file
            fprintf(outputs[part_index], "%s", buffer);
        }
    }

    for (i = 0; i < NUM_PARTS; i++) {
        fclose(outputs[i]);
    }
    fclose(fasta);
    free(buffer);

    if (verbose >= 1) {
        fprintf(stderr, "INFO : Done.\n");
    }

    return 0;
}
