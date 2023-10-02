///
/// @file
/// @brief Helper code for parsing JSON files used in MPCD inputs.
///
/// Variety of helper methods set up to aid in parsing JSON files for MPCD. Based around the cJson library released by
/// Dave Gramble, released under MIT license. Broadly speaking, these methods fall under one of the following
/// categories:
/// - Methods for handling the linked list containing all JSON tags/names recognised by the code, used for error
/// checking and verifying if any unknown tags have been submitted.
/// - Getting objects and values from the JSON file itself (i.e. ints, strings, arrays, etc.).
/// - File handling of the actual JSON input file.
///

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../headers/cJson.h"

// list of exclusions for "comment" tags
/// @brief The number of comment tags in the `commentTags` array.
const int commentTagCount = 4;
/// @brief The array of all possible comment tags.
char* commentTags[] = {"c", "comment", "//", "#"};

/* 
   Helper methods and structs to check if an element in the .json exists to the code or not
*/

///
/// @brief Initialise the linked list used for storing all found JSON tags.
///
/// Allocates memory required for the linked list, and sets the `next` and `str` members of the `head` to `NULL` and `""`.
///
/// @param head A pointer to the head of the linked list. Expected to be a `NULL` pointer passed as `&head`.
///
void initLL(linkedList **head){
   *head = (linkedList*) calloc(1, sizeof(linkedList));
   (*head)->next = NULL;
   dynAllocStr("", &((*head)->str)); // fill w blank string
}

///
/// @brief Debugging method to print the contents of the linked list to terminal.
///
/// Iterates through the linked list, printing each the `str` member of element to the terminal.
///
/// @param head The head of the linked list.
///
void printLL(linkedList *head){
   if (head == NULL) return;
   linkedList *curr = head;
   while(curr != NULL){
      printf("%s\n", curr->str);
      curr = curr->next;
   }
   
}

///
/// @brief Checks if a given string is present in the linked list. Returns 1 if true, 0 if false.
///
/// Iterates through the linked list and compares the `str` member of each element to the given string, returning 1 if
/// a match is found. If it is not found in the list, returns 0 at the end.
///
/// @param head The head of the linked list.
/// @param val The string value to check for.
/// @return Returns 1 if `val` is present in linked list, 0 if not.
///
int compareList(linkedList *head, const char* val){
  linkedList *current = head;
  while (current != NULL){ // iterate through linked list until end
     if (!strcmp(current->str, val)){
        return 1; //if we find something w the correct val then return true
     }

     current = current->next; //iterate
  }

  return 0; // if you're here then it doesnt exist in the list
}

///
/// @brief Checks if a given string is considered a valid "comment" tag. Returns 1 if true, 0 if false.
///
/// Iterates through the array of comment tags and compares the given string to each element, returning 1 if found. If
/// not found, it returns 0.
///
/// @param tag The string to check if it is a comment or not.
/// @return 1 if `tag` is a comment tag, 0 if not.
int isComment(const char* tag){
   int i; // counter
   for (i = 0; i < commentTagCount; i++) {
      if (strcmp(tag, commentTags[i]) == 0) {
          return 1;
      }
   }

   return 0; // if you get here then the tag isn't in the comment list
}

///
/// @brief Pushes a new element to the end of the linked list.
///
/// Iterates through the linked list until the end is found (ie, the `next` member of the current element is `NULL`).
/// Then allocates memory for a new element, and sets the `next` member of the last linked list node to point to it. The
/// `str` member of the new element is set to the given string value, and the `next` member of the new element is set to
/// NULL.
///
/// @param head The head of the linked list.
/// @param val The string value to be added as a new element to the end of the linked list.
///
void pushLL(linkedList * head, const char* val){
   linkedList *curr = head;
   while (curr->next != NULL){
      curr = curr->next; //iterate to end of list
   }

   // alloc and add to end
   curr->next = (linkedList*) calloc(1, sizeof(linkedList));
   dynAllocStr(val, &curr->next->str);
   curr->next->next = NULL;
}

///
/// @brief Frees the memory allocated to the linked list.
///
/// Iterates through the linked list, freeing the memory allocated to each element iteratively until the end of the list
/// is reached. The memory allocated to the `str` member of each element is freed first, then the memory allocated to
/// the linked list node itself.
///
/// @param head The head of the linked list.
///
void freeLL(linkedList *head){
   if (head == NULL) return;

   linkedList *curr = head;
   do
   {
      linkedList *next = curr->next; //save for later
      free(curr->str);
      free(curr);
      curr = next;
   } while (curr != NULL); // specifically use a do-while loop here to free all

   return;
}

///
/// @brief Verifies if there are any elements in the JSON file that are not recognised/read by MPCD.
///
/// A method to verify if there any "unknown" JSON tags in the given jObj. Method is recursive for arrays. It is assumed
/// that all parsing is done before this method is called (as parsing methods populate the `jsonTagList` linked list).
/// Calls `exit(EXIT_FAILURE)` if an unknown tag is found.
///
/// @param jObj The cJSON json object to check.
/// @param jsonTagList A linked list containing all JSON tags/ names recognised by MPCD.
/// @param arrayList A list of all array contents in the JSON file.
///
void verifyJson(cJSON *jObj, linkedList *jsonTagList, linkedList* arrayList){
   cJSON *childObj = NULL;
   int validJson = 1; // flag on whether the JSON is valid or not.

   cJSON_ArrayForEach(childObj, jObj){ // loop through children
      const char *jTag = childObj->string; // get json tag
      
      // check if this tag exists on the list of known objects
      if (!compareList(jsonTagList, jTag)){
         if (!isComment(jTag)){ // check if the tag is a comment
            // throw a warning if it's not
            printf("JSON Read Warning: Found unrecognised json tag: %s. Tag will be ignored.\n", jTag);
            validJson = 0;
         }
      } else {
         // if tag exists, check if it is an array and verify the subarray if necessary
         if (compareList(arrayList, jTag)) {

            // assume the array is of a child objects, which need to 
            //    be manually parsed
            cJSON *arrayChild = NULL;
            cJSON_ArrayForEach(arrayChild, childObj){
               verifyJson(arrayChild, jsonTagList, arrayList);
            }
         }
      }
   }

   if (!validJson){
      printf("JSON Read Error: Errors found.\n");
      exit(EXIT_FAILURE);
   }
}

/* 
   A variety of helper objects to read from Json and fill up a given data structure, and have some 
      mild segfault protection.
   These can be used to access elements within any given json object (ie, something surrounded by 
      {} )
   These SHOULD also be used to parse data from within array contained objects.
*/

///
/// @brief Set up a cJSON array object ready for parsing.
///
/// Searches for an array given by tag `val` in the `jObj` cJson object, and prepares a new cJson object `toReturn`
/// containing the child elements of the array. Keeps track of the `val` name for the json tag linked list, and also
/// marks the tag as one corresponding to an array of objects by pushing to the linked list `arrayList`.
///
/// This is set up to handle both primitives and objects in the array. Set `type` to 0 for primitives, and 1 for custom
/// objects.
///
/// @param jObj The cJson object to parse.
/// @param toReturn A pointer to the cJSON array object to be set up.
/// @param val The JSON tag of the array to be parsed.
/// @param jsonList A linked list containing all JSON tags/ names recognised by MPCD.
/// @param arrayList A linked list of all array objects in the JSON file.
/// @param type The type of array to be parsed. 0 for primitives (int, double, str, etc), 1 for custom objects (species, BCs, etc).
///
void getCJsonArray(cJSON *jObj, cJSON **toReturn, const char* val, 
        linkedList *jsonList, linkedList *arrayList, int type){
   *toReturn = cJSON_GetObjectItemCaseSensitive(jObj, val);
   pushLL(jsonList, val);
   if (type)   pushLL(arrayList, val); // if array of objects add to array list
}

///
/// @brief Parsing method for reading a primitive int from a cJSON object.
///
/// Returns an integer object from the given cJson file searching for a particular jsonTag. If no appropriate JSON tag
/// is found then it will return default value `d`.
///
/// @param cJSONRoot The root cJSON object to parse.
/// @param jsonTag The JSON tag/ name to search for in the object. Case-sensitive.
/// @param d The default value to return if the tag is not found.
/// @param head The head of the linked list containing all JSON tags/ names recognised by MPCD.
/// @return The parsed value of the int object corresponding to `jsonTag`.
///
int getJObjInt(cJSON *cJSONRoot, const char* jsonTag, int d, linkedList *head){
   pushLL(head, jsonTag); // add jsonTag to head

   // cJson bits
   cJSON *jObj = cJSON_GetObjectItemCaseSensitive(cJSONRoot, jsonTag);
   if (jObj == NULL){ // check for non-existence
      return d; 
   }
   
   int buff = jObj->valueint; // make a buffer to return an appropriate val
   return buff;
}

///
/// @brief Parsing method for reading a primitive int that can have multiple corresponding tags from a cJSON object.
///
/// Iterates through a list of potential JSON tags for an int object, checking if it is present using getJObjInt(). If
/// no appropriate JSON tag is found then it will return default value `d`.
///
/// @param cJSONRoot The root cJSON object to parse.
/// @param jsonTags An array of JSON tags/ names to search for in the object. Case-sensitive. The last tag is prioritised.
/// @param count The number of possible names for the JSON tag, ie the length of the `jsonTags` array.
/// @param d The default value to return if the tag is not found.
/// @param head The head of the linked list containing all JSON tags/ names recognised by MPCD.
/// @see getJObjInt()
/// @return The parsed value of the int object corresponding to one of `jsonTags`.
///
int getJObjIntMultiple(cJSON *cJSONRoot, const char** jsonTags, int count, int d, linkedList *head) {
    int i;
    int buff = d; // buffer to return appropriate val, set to default initially

    for (i = 0; i < count; i++) { // loop through tags
        const char* jsonTag = jsonTags[i]; // get current tag we iterate over

        int tagValue = getJObjInt(cJSONRoot, jsonTag, d, head); // get value of tag
        if (tagValue != d) { // if it's not the default value then we found something
            buff = tagValue; // set buffer to this value
        }
    }

    return buff;
}

///
/// @brief Parsing method for reading a primitive double from a cJSON object.
///
/// Returns a double object from the given cJson file searching for a particular jsonTag. If no appropriate JSON tag is
/// found then it will return default value `d`.
///
/// @param cJSONRoot The root cJSON object to parse.
/// @param jsonTag The JSON tag/ name to search for in the object. Case-sensitive.
/// @param d The default value to return if the tag is not found.
/// @param head The head of the linked list containing all JSON tags/ names recognised by MPCD.
/// @return The parsed value of the double object corresponding to `jsonTag`.
///
double getJObjDou(cJSON *cJSONRoot, const char* jsonTag, double d, linkedList *head){
   pushLL(head, jsonTag); // add jsonTag to head

   // cJson bits
   cJSON *jObj = cJSON_GetObjectItemCaseSensitive(cJSONRoot, jsonTag);
   if (jObj == NULL){ // check for non-existence
      return d; 
   }
   
   double buff = jObj->valuedouble; // make a buffer to return an appropriate val
   return buff;
}

///
/// @brief Helper function to dynamically allocate a string with value val to `toReturn`.
///
/// Allocates memory and copies the string `val` to `toReturn`.
///
/// @param val The value to be copied to the string.
/// @param toReturn The string object to be allocated to. Expects an object of form `&myStr`.
///
void dynAllocStr(const char *val, char **toReturn){
   int len = strlen(val) + 1; // length of memory to allocate
   *toReturn = malloc( len);
   if (*toReturn == NULL){ // stupid error checking
      printf("Failed to allocate memory for string: %s \n", val);
      return;
   }

   *toReturn = strcpy(*toReturn, val);
}

///
/// @brief Parsing method for reading a primitive string from a cJSON object.
///
/// Returns an string object from the given cJson file searching for a particular jsonTag. If no appropriate JSON tag
/// is found then it will return default value `d`. Unlike other methods, this is not leak-safe, and the user must free
/// it themselves.
///
/// @param cJSONRoot The root cJSON object to parse.
/// @param jsonTag The JSON tag/ name to search for in the object. Case-sensitive.
/// @param d The default value to return if the tag is not found.
/// @param toReturn String pointer, acts as a return pointer. Returns the string corresponding to the JSON tag specified.
/// @param head The head of the linked list containing all JSON tags/ names recognised by MPCD.
///
void getJObjStr(cJSON *cJSONRoot, const char* jsonTag, const char* d, char **toReturn, linkedList *head){
   pushLL(head, jsonTag); // add jsonTag to head

   // cJson bits
   cJSON *jObj = cJSON_GetObjectItemCaseSensitive(cJSONRoot, jsonTag);
   if (jObj == NULL){ // check for non-existence
      dynAllocStr(d, toReturn); // allocate the string dynamically
      return; 
   }
   
   const char* buff = jObj->valuestring; // make a buffer to return an appropriate val
   dynAllocStr(buff, toReturn); // allocate the string dynamically
   return;
}

///
/// @brief Gets a full string of the contents of a file to be parsed into a cJSON object afterwards.
///
/// Opens a file and reads the contents into a single string. The string is dynamically allocated to match the size of
/// the text in the file. If successful, returns 0, otherwise if there are errors doing this process an error code is
/// provided:
/// - 1: Error opening the file
/// - 2: Error allocating memory for the string
/// - 3: Error reading the file
///x
/// @param inFile The file pointer to be read from.
/// @param fileStr A string that will contain the contents of the file. Expecting to be passed an object of form
/// `&myStr`.
/// @return Returns 0 if successful, otherwise returns a pseudo-error code.
///
int getFileStr(char* inFile, char** fileStr){

	// printf("Reading file %s \n", inFile);
   
   FILE *fptr;
   // read file in with basic error checking
   if ((fptr = fopen(inFile, "r")) == NULL){ 
      printf("Error opening file %s\nQuitting.\n", inFile);
      return 1;
   }

   /* Get the number of bytes */
   fseek(fptr, 0L, SEEK_END);
   long numbytes = ftell(fptr) + 1;
   
   /* reset the file position indicator to 
   the beginning of the file */
   fseek(fptr, 0L, SEEK_SET);	
   
   /* grab sufficient memory for the 
   buffer to hold the text */
   *fileStr = (char*)calloc(numbytes, sizeof(char));	
   
   /* memory error */
   if(*fileStr == NULL){
      printf("Error allocating memory for file %s\n", inFile);
      return 2;
   }  
   
   /* copy all the text into the buffer */
   if (fread(*fileStr, sizeof(char), numbytes, fptr) == 0){
      printf("Could not read file %s\n", inFile);
      return 3;
   }

   fclose(fptr); // close file to free memory
   return 0;
}