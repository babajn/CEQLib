//#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <tgmath.h>
#include <assert.h>
//#include <mpi.h>
#include <sys/stat.h>
#include "Utilities.h"

bool CheckFile(const char FileName[]){
    FILE *fp;
    fp = fopen(FileName,"r");
    if(fp == NULL){
        return false;
    }
    fclose(fp);
    return true;
}

bool CheckDir(const char DirName[]){
    FILE *fp;
    fp = fopen(DirName,"r");
    if(fp == NULL){
        return false;
    }
    fclose(fp);
    return true;
}

void MakeDir(const char DirName[]){

    //if(MPIGetMyID() == MPI_ROOT_RANK){
        if(CheckDir(DirName) == false){
            int checkflag = mkdir(DirName,0755);
            if(checkflag < 0){
                fprintf(stderr,"Directory [%s] creation error.\n",DirName);
                exit(0);
            } else {
                fprintf(stderr,"Directory [%s] is created.\n",DirName);
            }
        } 
    //}

    return ;
}
