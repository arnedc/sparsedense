#include <stdio.h>
#include <stdlib.h>
#include <iosfwd>
#include <string>
#include "shared_var.h"
#include <iostream>
#include <fstream>
using namespace std;

int read_input ( char * filename ) {
    std::ifstream inputfile ( filename );
    string line;
    bool obs_bool=false, dense_bool=false, rand_bool=false, fixed_bool=false, Xfile_bool=false, Tfile_bool=false, Zfile_bool=false, lam_bool=false;
    bool blocksize_bool=false;

    filenameT= ( char* ) calloc ( 100,sizeof ( char ) );
    filenameX= ( char* ) calloc ( 100,sizeof ( char ) );
    filenameZ= ( char* ) calloc ( 100,sizeof ( char ) );
    filenameD= ( char* ) calloc ( 100,sizeof ( char ) );
    filenameC= ( char* ) calloc ( 100,sizeof ( char ) );

    lambda=100;
    blocksize=64;
    printD_bool=false;
    printsparseC_bool=false;


    while ( std::getline ( inputfile,line ) ) {
        if ( line=="#Observations" ) {
            std::getline ( inputfile,line );
            n=atoi ( line.c_str() );
            obs_bool=true;
        } else if ( line=="#DenseEffects" ) {
            std::getline ( inputfile,line );
            k=atoi ( line.c_str() );
            dense_bool=true;
        } else if ( line=="#FixedEffects" ) {
            std::getline ( inputfile,line );
            m=atoi ( line.c_str() );
            fixed_bool=true;
        } else if ( line=="#RandomEffects" ) {
            std::getline ( inputfile,line );
            l=atoi ( line.c_str() );
            rand_bool=true;
        } else if ( line=="#FileFixedEffects" ) {
            std::getline ( inputfile,line );
            line.copy ( filenameX,100 );
            Xfile_bool=true;
        } else if ( line=="#DenseDataFile" ) {
            std::getline ( inputfile,line );
            line.copy ( filenameT,100 );
            Tfile_bool=true;
        } else if ( line=="#RandDataFile" ) {
            std::getline ( inputfile,line );
            line.copy ( filenameZ,100 );
            Zfile_bool=true;
        } else if ( line=="#Lambda" ) {
            std::getline ( inputfile,line );
            lambda=atof ( line.c_str() );
            lam_bool=true;
        } else if ( line=="#BlockSize" ) {
            std::getline ( inputfile,line );
            blocksize=atoi ( line.c_str() );
            blocksize_bool=true;
	} else if ( line=="#OutputFileD" ) {
            std::getline ( inputfile,line );
            line.copy ( filenameD,100 );;
            printD_bool=true;
	} else if ( line=="#OutputFileSparseC" ) {
            std::getline ( inputfile,line );
            line.copy ( filenameC,100 );;
            printsparseC_bool=true;
        } else if ( line[0]=='/' || line.size() ==0 ) {}
        else {
            printf ( "Unknown parameter in inputfile, the following line was ignored: \n" );
            printf ( "%s\n",line.c_str() );
        }
    }
    if ( obs_bool ) {
        if ( dense_bool ) {
            if ( fixed_bool ) {
                if ( Tfile_bool ) {
                    if ( Xfile_bool ) {
                        if ( rand_bool ) {
                            if ( Zfile_bool ) {
                                if ( *position==0 && * ( position+1 ) ==0 ) {
                                    printf ( "number of observations:   \t %d\n", n );
                                    printf ( "number of dense effects:  \t %d\n", k );
                                    printf ( "number of fixed effects:  \t %d\n", m );
                                    printf ( "number of random effects: \t %d\n", l );
                                    printf ( "filename of dense effects:\t %s\n", filenameT );
                                    printf ( "filename of fixed effects:\t %s\n", filenameX );
                                    printf ( "filename of random effects:\t %s\n", filenameZ );
                                }
                            } else {
                                printf ( "ERROR: filename of random effects was not in input file or not read correctly\n" );
                                return -1;
                            }
                        } else {
                            printf ( "ERROR: number of random effects was not in input file or not read correctly\n" );
                            return -1;
                        }
                    } else {
                        printf ( "ERROR: filename of fixed effects was not in input file or not read correctly\n" );
                        return -1;
                    }
                } else {
                    printf ( "ERROR: filename of dense effects was not in input file or not read correctly\n" );
                    return -1;
                }
            } else {
                printf ( "ERROR: number of fixed effects was not in input file or not read correctly\n" );
                return -1;
            }
        } else {
            printf ( "ERROR: number of dense effects was not in input file or not read correctly\n" );
            return -1;
        }
    } else {
        printf ( "ERROR: number of observations was not in input file or not read correctly\n" );
        return -1;
    }
    if ( *position==0 && * ( position+1 ) ==0 ) {
        if ( blocksize_bool ) {
            printf ( "Blocksize of %d was used to distribute matrices across processes\n", blocksize );
        } else {
            printf ( "Default blocksize of %d was used to distribute matrices across processes\n", blocksize );
        }
        if ( lam_bool )
            printf ( "Start value of %g was used to estimate variance component lambda\n", lambda );
        else
            printf ( "Default start value of %g was used to estimate variance component lambda\n", lambda );
	if ( printD_bool )
	  printf("Matrix D will be written to binary file %s \n", filenameD);
	if ( printsparseC_bool )
	  printf("Sparse matrix C will be written in CSR format to text file %s \n", filenameC);
    }
    return 0;
}
