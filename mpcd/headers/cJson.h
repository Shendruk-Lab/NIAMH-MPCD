/*
   Helper code to set up JSON parsers in C.
   Uses the cJson library by Dave Gramble released under MIT license.
   Started by Timofey Kozhukhov, Summer 2021.
*/

#ifndef HGUARD_cJson_h
#define HGUARD_cJson_h

#include "../../dependencies/cJson/cJSON.h"

// a dumb linked list for handling list of saved strings
typedef struct linkedList linkedList;
struct linkedList{
   linkedList *next; 
   char* str;
};

/*
    Methods for handling the above linked list. Specifically designed for use 
        with the verifyJson method.
*/

void initLL(linkedList **head);
void printLL(linkedList *head);
int compareList(linkedList *head, const char* val);
void pushLL(linkedList * head, const char* val);
void freeLL(linkedList *head);

/*
    Methods for parsing and handling cJson data
*/

void getCJsonArray(cJSON *jObj, cJSON **toReturn, const char* val, 
   linkedList *jsonList, linkedList *arrayList, int type);
int getJObjInt(cJSON *cJSONRoot, const char* jsonTag, int d, linkedList *head);
int getJObjIntMultiple(cJSON *cJSONRoot, const char** jsonTags, int count, int d, linkedList *head);
double getJObjDou(cJSON *cJSONRoot, const char* jsonTag, double d, 
    linkedList *head);
void getJObjStr(cJSON *cJSONRoot, const char* jsonTag, const char* d, 
    char **toReturn, linkedList *head);
void verifyJson(cJSON *jObj, linkedList *jsonTagList, linkedList* arrayList);

// Helper methods for the above routines
void dynAllocStr(const char *val, char **toReturn);
int getFileStr(char* inFile, char** fileStr);

#endif