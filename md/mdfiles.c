
//================================================================================
//
// name:   mdfiles.c
// author: ftessier
// date:   2005-05-03 @ 11:04:36
//
//================================================================================
///
/// @file mdfile.c
/// Functions related to reading and writing of the data. Also contains a function
/// to parse the command-line arguments and some utility functions needed exclusively
/// for file operations.
///
//================================================================================


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <pwd.h>
#include "mdtypes.h"
#include "mderror.h"
#include "mdmemory.h"
#include "mdhistogram.h"
#include "mdsetup.h"
#include "mdsrd.h"
#include "mdfiles.h"
#include "../mpcd/headers/globals.h"


/// Parse the command-line options of the program using getopt. Options can be
/// saved in an options structure for later used in the parent function.
///
/// @param		argc argument count passed from main
/// @param		argv the argument values array
/// @param		sim a pointer to a simulation structure
/// @param		options a pointer to an options structure
/// @warning 	Don't use the more sophisticated argp functions because they are
/// 			not standard on other compilers; On the other hand, getopt
/// 			is part of POSIX.

//================================================================================
void ParseOptions (int argc, char *argv[], simptr sim, simoptions *options)
//================================================================================
{
	int c;
	int	opt_o=0, opt_i=0, opt_c=0;
	int arg;

	// get program name and pid
	strncpy (sim->programName, argv[0], STRMAX);
	sim->pid = getpid();

	// parse command-line options
	while ((c = getopt (argc, argv, "i:o:c:h:I:O:L")) != -1) {
		switch (c) {
			case 'i': // read input from the global input file string
				snprintf (sim->inputFile, STRMAX, "%s", mdInputFile);
				opt_i = 1;
				break;

			//FIXME: this is shit, doesnt take into account any further output. 
			case 'L': // legacy input option
				// name of input file
				snprintf (sim->inputFile, STRMAX, "%s/md.inp", optarg);
				// if (strcmp(sim->inputFile, optarg)) error (EPARSE);
				opt_i = 1;
				break;

			case 'o':
				// name of output directory
				snprintf (sim->outputDir, STRMAX, "%s", optarg);
				if (strcmp(sim->outputDir, optarg)) error (EPARSE);
				opt_o = 1;
				break;

			case 'c':
				// continue from checkpoint
				snprintf (sim->outputDir, STRMAX, "%s", optarg);
				if (strcmp(sim->outputDir, optarg)) error (EPARSE);
				opt_c = 1;
				options->setupType = SIMOPT_SETUP_CHKPOINT;
				break;

			case 'h':
				// help
				fprintf (stderr, "Usage:	md [options] -o <output directory>.\n");
				fprintf (stderr, "Options:  -i   name of input file, or stdin assumed\n");
				fprintf (stderr, "		  -c   name of checkpoint directory\n");
				exit(0);
				break;

			case ':':
				// missing argument
				fprintf (stderr, "Option `-%c' requires an argument.\n", optopt);
				error (EPARSE);
				break;

			case '?':
				// unknown option
				if (isprint (optopt))
				fprintf (stderr, "Unknown option `-%c'.\n", optopt);
				else
				fprintf (stderr, "Unknown option character `0x%x'.\n", optopt);
				// error (EPARSE);
				break;
			// default:
			// 	error (EPARSE);
			// 	break;
		}
	}

	// parse extra multichara cmd line options
	for(arg=1; arg<argc; arg++) {
		// Check for a dash
		if(argv[arg][0]=='-') {
			switch (argv[arg][1]) {
				case 'L': // legacy 
					if (argv[arg][2] == 'i'){ // arg 'Li' for legacy input
						arg++;
						
						snprintf (sim->inputFile, STRMAX, "%s/md.inp", argv[arg]);
						// if (strcmp(sim->inputFile, optarg)) error (EPARSE);
						opt_i = 1;
						break;
					}
			}
		}
	}

	// use stdin if no input file is specified
	if (!opt_c && !opt_i) {
		snprintf (sim->inputFile, STRMAX, "stdin");
	}

	// name of output directory is mandatory if sim not run from checkpoint
	if (!opt_c && !opt_o) {
		fprintf (stderr, "You need to specify an output directory with '-o'.\n");
		error (EPARSE);
	}
}


//================================================================================
void SetupDataFiles (simptr sim)
//================================================================================
{
	// Sets up the output files for the md simulation and opens them for
	// writing. Follow the template to add new files. Note that the log file
	// MUST have "log" as its label, while the checkpoint file MUST use
	// "chk". These could be specified by defines if it needs to be
	// flexible. For now it is hardcoded.

	fileptr	fptr;
	FILE	*pidfile;

	// log file
	fptr = FileNew (sim->files);
	snprintf (fptr->label,   	STRMAX, "log");
	snprintf (fptr->name,		STRMAX, "%s.log", sim->label);
	snprintf (fptr->desc,		STRMAX, "MD simulation message log");
	snprintf (fptr->type,		STRMAX, "ascii");
	snprintf (fptr->layout, 	STRMAX, "columns");
	snprintf (fptr->columns, 	STRMAX, "date, time (UTC), step, message");

	// record list in simulation structure
	sim->files = fptr;

	// checkpoint file
	fptr = FileNew (sim->files);
	snprintf (fptr->label,  	STRMAX, "chk");
	snprintf (fptr->name,		STRMAX, "%s.chk", sim->label);
	snprintf (fptr->desc,		STRMAX, "Checkpoint file");
	snprintf (fptr->type,		STRMAX, "binary");
	snprintf (fptr->layout,  	STRMAX, "blocks");
	snprintf (fptr->columns, 	STRMAX, "none");

	// include additional user file definitions
	#include "mdfilesdef.h"

	// open all data files
	fptr = sim->files;
	while (fptr) {

		// skip checkpoint file
		if (!strcmp(fptr->label, "chk")) {
			fptr = fptr->next;
			continue;
		}

		// open file for writing
		if (!(fptr->stream = fopen (fptr->name, "w"))) error (EFILE);

		// write data header at top of file
		WriteHeader (fptr, sim);

		// report in log
		LOG ("Created new simulation file %s\n", fptr->name);

		// next file
		fptr = fptr->next;
	}

	// write simulation process id (pid) in a file
	if (!(pidfile = fopen("pid", "w"))) error (EFILE);
	fprintf (pidfile, "%05d", sim->pid);
	fclose (pidfile);
}


//================================================================================
void ReopenDataFiles (simptr sim)
//================================================================================
{
	// Reopens data files when a chekpointed simulation resumes. The list of
	// files is stored in the checkpoint file and is loaded before this function
	// is called.

	fileptr	fptr;
	FILE	*pidfile;
	int read;

	// reopen files
	fptr = sim->files;
	while (fptr) {

		// skip checkpoint file
		if (!strcmp(fptr->label, "chk")) {
			fptr = fptr->next;
			continue;
		}

		// open file for appending
		if (!(fptr->stream = fopen (fptr->name, "a"))) error (EFILE);

		// truncate file to last checkpoint position
		fseek	  (fptr->stream, fptr->chkpos, SEEK_SET);
		read = ftruncate (fileno(fptr->stream), fptr->chkpos);
		if(read!=0) printf("Warning: MD Re-open read failed\n");

		// if we are reopening the log file
		if (!strcmp(fptr->label, "log")) LOG ("Restarting simulation from checkpoint\n");

		// report
		LOG ("Reopenning data file %s (offset %ld)\n", fptr->name, fptr->chkpos);

		// next file
		fptr = fptr->next;
	}

	// write simulation pid in a file
	if (!(pidfile = fopen("pid", "w"))) error (EFILE);
	fprintf (pidfile, "%05d", sim->pid);
	fclose (pidfile);
}


//================================================================================
void ReadParameters (char *inputFile, paramptr param, int nParam, char *label, char *stop)
//================================================================================
{
	// Reads nParam simulation parameters from the input file, enforcing the
	// correct parameter count. Parameters to be read are pointed by elements of
	// the paramList pointer list. Spaces and the order of the lines are
	// irrelevant in the input file. This function does not feature strong error
	// checking (i.e., it is easy to crash it with bad input format!).

	int			i, k, n, ok;
	char 		**lines, **strptr, *str;
	char		key[KEYSIZE], digits[]="eX0123456789.-\0";
	paramptr	ptr;
	FILE		*input;


	// open input file
	if (!strcmp(inputFile, "stdin")) {
		input = stdin;
	}
	else {
		input = fopen (inputFile, "r");
		if (!input) error (EINPUT);
	}

	// dump input lines into array
	lines = LinesToArray (input, stop);

	// close input file
	if (input != stdin) fclose (input);

	// filter unwanted characters and translate macros
	i=0;
	strptr = &lines[i];
	while (*strptr) {
		FilterString (strptr, " \t");
		TranslateMacros (strptr);
		strptr = &lines[++i];
	}

	// loop over all parameters
	for (n=0; n<nParam; n++) {

		// get pointer to parameter
		ptr = param+n;

		// construct the key string for the current parameter
		snprintf (key, KEYSIZE, "%s%s=", label, ptr->name);

		// search for the key in all the input lines, up to stop line
		i=0; ptr->n=0;
		while (!strstr (str=lines[i++], stop)) {

			// check if this line begins with the sought key
			if (!strncmp (str, key, strlen(key))) {

			// skip to after equal sign
			str = str+strlen(key);

			// read the value(s) for the parameter
			for (k=0; k<ptr->count; k++) {

				// read value(s) depending on type
				switch (ptr->type) {
					case INTG:
						str = strpbrk (str, digits);
						if (str) ptr->n += sscanf (str, "%d",  (int *) ptr->value + k);
						break;
					case HEXA:
						str = strpbrk (str, digits);
						if (str) ptr->n += sscanf (str, "%X",  (int *) ptr->value + k);
						break;
					case REAL:
						str = strpbrk (str, digits);
						if (str) ptr->n += sscanf (str, real_FORMAT_STR, (real *) ptr->value + k);
						break;
				}

				// advance to next item
				if (str) str = str + strspn (str, digits);
				else break;
			}
			break;
			}
		}
	}

	// verify parameters
	ok = 1;
	for (n=0; n<nParam; n++) {
		ptr = param+n;
		if (ptr->n != ptr->count) {
			fprintf (stderr, "\nError reading parameter: %s\n", ptr->name);
			fprintf (stderr, "%d value(s) required, but only %d found!\n", ptr->count, ptr->n);
			ok = 0;
		}
	}
	if (!ok) error (EABORT);

	// free memory used to store the input lines
	i=0;
	while (lines[i]) free (lines[i++]);
	free (lines);
}


//================================================================================
void WriteParameters (FILE *output, paramptr param, int nParam, char *label)
//================================================================================
{
	// Write a list of nParam simulation parameters to the output
	// file. Parameters to be written are in the passed param list.

	int			k, n;
	paramptr	ptr;

	// write parameter list to file
	for (n=0; n<nParam; n++) {
		ptr = param+n;

		// parameter name
		fprintf (output, "%s             %-20s = ", label, ptr->name);

		// openning paranthesis for vector parameters
		if (ptr->count > 1) fprintf (output, "(");

		// parameter value(s)
		for (k=0; k < ptr->count; k++) {

			// write value depending on type
			switch (ptr->type) {
				case INTG:
					fprintf (output, "%d", *(((int*)ptr->value)+k));
					break;
				case HEXA:
					if (*(((unsigned int*)ptr->value) + k) == 0) fprintf (output, "0");
					else fprintf (output, "%#010X", *(((unsigned int*) ptr->value) + k));
					break;
				case REAL:
					fprintf (output, "%G", *(((real*) ptr->value) + k));
					break;
			}

			// comma for vector parameters
			if (k >= 0 && k < ptr->count-1) fprintf (output, ", ");
		}

		// closing paranthesis for vector parameters
		if (ptr->count > 1) fprintf (output, ")");

		// trailing newline
		fprintf (output, "\n");
	}

	fprintf (output, "# \n");
}


//================================================================================
char **LinesToArray (FILE *input, const char *stop)
//================================================================================
{
	// Transfers the lines of an input stream to an array of strings, up to the
	// line containing the stop string. Returns a pointer to the NULL
	// terminated array of strings. This function works for an arbitrary number
	// of lines of arbitrary lengths (limited only by available memory).

	char **lines;
	int  k, numlines, numchars, done;

	// allocate a block of lines
	lines = (char **) mycalloc (LINEBLOCK+1, sizeof(char *));
	numlines = LINEBLOCK;

	k=0;
	do {
		// allocate memory for additional line pointers if needed
		if (k >= numlines) {
			numlines += LINEBLOCK;
			lines = (char **) myrealloc (lines, (numlines+1)*(sizeof(char *)));
		}

		// make sure the lines array is always NULL terminated
		lines[k+1] = NULL;

		// allocate memory for a line
		lines[k] = (char *) mycalloc (LINESIZE+1, sizeof(char));

		// read one line of arbitrary length
		numchars = LINESIZE;
		do {
			fgets (lines[k]+numchars-LINESIZE, LINESIZE+1, input);
			done = (strchr(lines[k],'\n') || feof(input));
			// allocate space for more characters if needed
			if (!done) {
				numchars += LINESIZE;
				lines[k] = (char *) myrealloc (lines[k], (numchars+1)*(sizeof(char)));
			}

		} while (!done);

		// stop if we reached the stop string
		if (strstr(lines[k], stop)) break;

		// increment line counter
		k++;

	} while (!feof(input));

	return lines;
}


//================================================================================
void FilterString (char **stringptr, const char *filter)
//================================================================================
{
	// Filters characters out of a string (e.g. spaces, tabs, etc.); the
	// original string is modified.

	char *ptr1=*stringptr, *ptr2;

	// remove unwanted characters from string
	while ((ptr1 = strpbrk (ptr1, filter))) {
		ptr2 = ptr1 + strspn (ptr1, filter);	// ptr2 -> next non-filtered char
		memmove (ptr1, ptr2, strlen(ptr2)+1);	// move *ptr2 to *ptr1
	}
}


//================================================================================
void TranslateMacros (char **stringptr)
//================================================================================
{
	// Translate macros found in the input string; the original string is
	// modified. Note that the string pointer may also be modified if the string
	// is not long enough to hold the translated string.

	char macro[STRLEN], value[STRLEN];

	macro[STRMAX]='\0';
	value[STRMAX]='\0';

	snprintf (macro, STRMAX, "TYPE_WALL");
	snprintf (value, STRMAX, "%u", TYPE_WALL);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "TYPE_FLUID");
	snprintf (value, STRMAX, "%u", TYPE_FLUID);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "TYPE_MONOMER");
	snprintf (value, STRMAX, "%u", TYPE_MONOMER);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "LAYOUT_FLUID");
	snprintf (value, STRMAX, "%u", LAYOUT_FLUID);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "LAYOUT_WALL");
	snprintf (value, STRMAX, "%u", LAYOUT_WALL);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "LAYOUT_SURFACE");
	snprintf (value, STRMAX, "%u", LAYOUT_SURFACE);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "LAYOUT_ANCHOR");
	snprintf (value, STRMAX, "%u", LAYOUT_ANCHOR);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "LAYOUT_PLATES");
	snprintf (value, STRMAX, "%u", LAYOUT_PLATES);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "LAYOUT_CYLINDER");
	snprintf (value, STRMAX, "%u", LAYOUT_CYLINDER);
	ReplaceMacro (stringptr, macro, value);

	// Tyler added the following
	snprintf (macro, STRMAX, "LAYOUT_RODX");
	snprintf (value, STRMAX, "%u", LAYOUT_RODX);
	ReplaceMacro (stringptr, macro, value);

	// Zahra added the following
	snprintf (macro, STRMAX, "LAYOUT_RODY");
	snprintf (value, STRMAX, "%u", LAYOUT_RODY);
	ReplaceMacro (stringptr, macro, value);

	// Karolina added the following
	snprintf (macro, STRMAX, "LAYOUT_U");
	snprintf (value, STRMAX, "%u", LAYOUT_U);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "GROUP_NONE");
	snprintf (value, STRMAX, "%#010X", GROUP_NONE);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "GROUP_FLUID");
	snprintf (value, STRMAX, "%#010X", GROUP_FLUID);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "GROUP_WALL_OUT");
	snprintf (value, STRMAX, "%#010X", GROUP_WALL_OUT);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "GROUP_WALL_IN");
	snprintf (value, STRMAX, "%#010X", GROUP_WALL_IN);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "GROUP_WALL");
	snprintf (value, STRMAX, "%#010X", GROUP_WALL);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "GROUP_ION");
	snprintf (value, STRMAX, "%#010X", GROUP_ION);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "GROUP_ION_POS");
	snprintf (value, STRMAX, "%#010X", GROUP_ION);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "GROUP_ION_NEG");
	snprintf (value, STRMAX, "%#010X", GROUP_ION);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "GROUP_MONOMER");
	snprintf (value, STRMAX, "%#010X", GROUP_ION);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "GROUP_ALL");
	snprintf (value, STRMAX, "%#010X", GROUP_ALL);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "GEOM_BULK");
	snprintf (value, STRMAX, "%u", GEOM_BULK);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "GEOM_CAPILLARY");
	snprintf (value, STRMAX, "%u", GEOM_CAPILLARY);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "GEOM_PLATES");
	snprintf (value, STRMAX, "%u", GEOM_PLATES);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "GEOM_CYLINDER");
	snprintf (value, STRMAX, "%u", GEOM_CYLINDER);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "LATT_SC");
	snprintf (value, STRMAX, "%u", LATT_SC);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "LATT_BCC");
	snprintf (value, STRMAX, "%u", LATT_BCC);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "LATT_FCC");
	snprintf (value, STRMAX, "%u", LATT_FCC);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "LATT_BULK");
	snprintf (value, STRMAX, "%u", LATT_BULK);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "CUTOFF_WCA");
	snprintf (value, STRMAX, "%d", CUTOFF_WCA);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "AUTO");
	snprintf (value, STRMAX, "%d", -1);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "SURFACE");
	snprintf (value, STRMAX, "%d", SURFACE);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "VOLUME");
	snprintf (value, STRMAX, "%d", VOLUME);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "NONE");
	snprintf (value, STRMAX, "%d", NONE);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "YES");
	snprintf (value, STRMAX, "%d", YES);
	ReplaceMacro (stringptr, macro, value);

	snprintf (macro, STRMAX, "NO");
	snprintf (value, STRMAX, "%d", NO);
	ReplaceMacro (stringptr, macro, value);
}


//================================================================================
void ReplaceMacro (char **stringptr, char *macro, char *value)
//================================================================================
{
	// Replaces macros in the line with their actual values. If the expanded
	// string is larger than the original one, we realloc the string pointer to
	// ensure we have enough space. Use a while loop to make sure we catch all
	// instance of the macro in each line.

	int	 length;
	char *ptr1, *ptr2;

	while ((ptr1 = strstr (*stringptr, macro))) {

		// predict length of translated string
		length = strlen(*stringptr) - strlen(macro) + strlen(value);

		// make sure we have enough space for the translation
		if (length > strlen(*stringptr)) {
			*stringptr = (char *) myrealloc (*stringptr, (length+1)*(sizeof(char)));
			ptr1 = strstr (*stringptr, macro);
		}

		// get the remaider of the line
		ptr2  = ptr1 + strlen(macro);

		// shift remainder of line to appropriate position
		memmove (ptr1+strlen(value), ptr2, strlen(ptr2)+1);

		// write translated macro value
		strncpy (ptr1, value, strlen(value));
	}
}


//================================================================================
void SetSimLabel (char *label, int pid)
//================================================================================
{
	// Sets the file label from the date and the pid of the simulation process.

	time_t		seconds;
	struct tm	*t;

	// get current time and convert to UTC
	time (&seconds);
	t = gmtime (&seconds);

	// set file label
	snprintf (label, STRMAX, "%04d%02d%02d-%05d", 1900+t->tm_year, 1+t->tm_mon, t->tm_mday, pid);
}


//================================================================================
void SetWorkingDir (char *dir, char *simdir)
//================================================================================
{
	// Changes the current working directory to that stored in sim->outputDir
	// and creates a new directory for the simulation output files.

	// go to specified directory
	if (chdir (dir) == -1) error (ERRNO);

	// make a directory for this simulation
	if (mkdir(simdir, 0700) == -1) error (ERRNO);

	// go to newly created directory
	if (chdir (simdir) == -1) error (ERRNO);
}


/// Print a short message with timestamp to the specified stream. This is useful
/// for logging information about the simulation. The function first prints a time
/// stamp and then the message passed in msg. This function is usually called
/// through the LOG variadic macro, which can handles a format string and a variable
/// number of arguements.
///
/// @param		sim a pointer to a simulation structure
/// @param		stream the file stream where the message should be printed
/// @param		msg the message to print
/// @return		void
/// @warning 	This function uses the overall step count in the sim structure to label
///             the entries, after the timestamp
/// @see		The LOG variadic macro

//================================================================================
void Report (simptr sim, FILE *stream, char *msg)
//================================================================================
{
        // Prints a timestamp and message to the specified stream

        time_t          seconds;
        struct tm       *t;

        // get current time and convert to UTC
        time (&seconds);
        t = gmtime (&seconds);

        // print time stamp and flag
        fprintf (stream, "%04d-%02d-%02d %02d:%02d:%02d %09d ",
				 1900+t->tm_year, 1+t->tm_mon, t->tm_mday,
				 t->tm_hour, t->tm_min, t->tm_sec, sim->step[count_]);

        // print message
	fprintf (stream, "%s",msg);
        fflush  (stream);
}


//================================================================================
void WriteHeader (simfile *f, simptr sim)
//================================================================================
{
	// Writes a standard header containing meta-information about the data file
	// and the simulation itself. This header should be included at the
	// beginning of every file pertaining to the simulation. In particular, all
	// the parameters passed to the simulation are written in this header.

	pid_t 			pid;
	time_t			seconds;
	char  			hostname[64], *username, na[]= "n/a";
	struct tm 		*t;
	struct passwd	*user;

	// get system and user information
	gethostname (hostname, 64);
	pid  = getpid();
	if ((user = getpwuid (geteuid())))
		username = user->pw_name;
	else
		username = na;

	// get current time and convert to UTC
	time (&seconds);
	t = gmtime (&seconds);

	// section: content
	fprintf (f->stream, "##############################################################################################\n");
	fprintf (f->stream, "#\n");
	fprintf (f->stream, "# section:               content\n");
	fprintf (f->stream, "# content-id:            %04d%02d%02d-%02d%02d%02d-%05d-%s-%s\n", \
				1900+t->tm_year, 1+t->tm_mon, t->tm_mday, t->tm_hour, t->tm_min, \
				t->tm_sec, (int) pid, hostname, f->label);
	fprintf (f->stream, "# content-owner:         %s\n", username);
	fprintf (f->stream, "# content-desc:          %s\n", f->desc);
	fprintf (f->stream, "# content-type:          %s\n", f->type);
	fprintf (f->stream, "# content-layout:        %s\n", f->layout);
	fprintf (f->stream, "# content-columns:       %s\n", f->columns);
	fprintf (f->stream, "# \n");

	// section: parameter
	fprintf (f->stream, "# section:               parameter\n");
	WriteParameters (f->stream, sim->param, sim->nParam, "# parameter:");

	// section: program
	fprintf (f->stream, "# section:               program\n");
	fprintf (f->stream, "# program-name:          %s\n", sim->programName);
	fprintf (f->stream, "# program-author:        Frederic Tessier\n");
	fprintf (f->stream, "# \n");

	// end-header
	fprintf (f->stream, "##############################################################################################\n");
	fprintf (f->stream, "# end-header\n");

	// flush stream
	fflush (f->stream);
}


//================================================================================
void CheckpointWrite (void *simvoid)
//================================================================================
{
	// Saves the current simulation state to disk so that it can be restarted
	// from this preserved state at a later time. Checkpoint files are written
	// in binary format for speed and simplicity. Be careful: different
	// platforms may have different byte ordering, and checkpoint files from
	// different program versions will not in general be compatible.

	char		tmpName[STRLEN];
	fileptr		fptr, chkptr;
	histptr		hptr;
	sceneptr	sptr;
	FILE		*chk;
	simptr		sim = (simptr) simvoid;

	// get the checkpoint simfile
	chkptr = GetSimFile (sim->files, "chk");
	if (!chkptr) error (ECHKWRITE);

	// copy previous checkpoint file to temporary file
	tmpName[STRMAX]='\0';
	snprintf (tmpName, STRMAX, "%s~", chkptr->name);
	rename   (chkptr->name, tmpName);

	// report
	LOG ("*** CHECKPOINT ***\n");

	// open new checkpoint file for writing.
	if (!(chkptr->stream = fopen (GetSimFileName(sim->files, "chk"), "w"))) error (EFILE);
	chk = chkptr->stream;

	// get checkpoint position of all files
	fptr = sim->files;
	while (fptr) {
		fptr->chkpos = ftell (fptr->stream);
		fptr = fptr->next;
	}

	// update neighbor list
	RefreshNebrList (sim);

	// write simulation data (sim strucutre must be first)
	WriteBinary (BINARY_DESC_SIM,			1,			 	sizeof(simulation),  	(void *) sim, chk);
	WriteBinary (BINARY_DESC_ATOM,			sim->atom.n,	sizeof(particleMD),		(void *) sim->atom.items, 	 chk);
	WriteBinary (BINARY_DESC_ATOM_START,	1,			 	sizeof(particleMD *),  	(void *) &(sim->atom.items), chk);
	WriteBinary (BINARY_DESC_POLYMER,		sim->polymer.n, sizeof(listPoly),		(void *) sim->polymer.items, chk);

	// write simfiles
	fptr = sim->files;
	while (fptr) {
		// simfile record
		WriteBinary (BINARY_DESC_FILE, 1, sizeof(simfile), (void *) fptr, chk);
		fptr = fptr->next;
	}

	// write histograms
	hptr = sim->histograms;
	while (hptr) {
		// histogram record
		WriteBinary (BINARY_DESC_HISTOGRAM, 1, sizeof(simhist), (void *) hptr, chk);
		// data
		WriteBinary (BINARY_DESC_HISTBIN,   hptr->n, sizeof(real), (void *) hptr->bin,   chk);
		WriteBinary (BINARY_DESC_HISTCOUNT, hptr->n, sizeof(real), (void *) hptr->count, chk);
		hptr = hptr->next;
	}

	// write scenes
	sptr = sim->scenes;
	while (sptr) {
		// scene record
		WriteBinary (BINARY_DESC_SCENE, 1, sizeof(simscene), (void *) sptr, chk);
		sptr = sptr->next;
	}

	// close checkpoint file
	fclose (chk);

	// remove temporary file
	unlink (tmpName);
}


//================================================================================
void CheckpointRead (simptr sim)
//================================================================================
{
	// Reads a checkpointed simulation state to restart a simulation. The
	// checkpoint file is in binary format, and each data block is preceded by
	// binary int description and count fields. We read each block according to
	// the these values.

	int			i, n, desc, count, len;
	char 		*basename;
	particleMD	*atomStart;
	fileptr		*fpptr, f;
	histptr		*hpptr, h;
	sceneptr	*spptr, s;
	FILE		*chk;

	// set up lists
	sim->files = NULL;
	sim->histograms = NULL;
	fpptr = &(sim->files);
	hpptr = &(sim->histograms);
	spptr = &(sim->scenes);

	// remove trailing slash(es) in directory name
	len = strlen (sim->outputDir);
	while (sim->outputDir[len-1] == '/') sim->outputDir[len-1] = '\0';

	// go inside output directory
	if (chdir (sim->outputDir) == -1) error (ERRNO);

	// get base name of file
	basename = strrchr (sim->outputDir, '/');
	if (!basename) basename = sim->outputDir;
	else basename += 1;

	// construct checkpoint file name
	snprintf (sim->inputFile, STRMAX, "%s.chk", basename);

	// open checkpoint file for reading
	if (!(chk = fopen (sim->inputFile, "r"))) error (ECHKREAD);

	// read binary checkpoint data
	i=0;
	while (!feof(chk)) {

		// read data description and count
		n = fread (&desc,  sizeof(int), 1, chk); if (n != 1 && !feof(chk)) error (ECHKREAD);
		n = fread (&count, sizeof(int), 1, chk); if (n != 1 && !feof(chk)) error (ECHKREAD);
		if (feof(chk)) break;

		// must read simulation structure first
		if (i==0 && desc != BINARY_DESC_SIM) error (ECHKREAD);
		i++;

		// read data depending on description
		switch (desc) {

			case (BINARY_DESC_SIM):
				// simulation structure
				n = fread (sim, sizeof(simulation), count, chk); if (n != count) error (ECHKREAD);
				break;

			case (BINARY_DESC_ATOM):
				// atoms
				sim->atom.items = (particleMD *) mycalloc (sim->atom.max, sizeof(particleMD));
				n = fread (sim->atom.items, sizeof(particleMD), count, chk); if (n != count) error (ECHKREAD);
				break;

			case (BINARY_DESC_ATOM_START):
				// saved atom base pointer
				n = fread (&atomStart, sizeof(particleMD *), count, chk); if (n != count) error (ECHKREAD);
				break;

			case (BINARY_DESC_POLYMER):
				// polymers
				sim->polymer.items = (itemPoly *) mycalloc (sim->polymer.max, sizeof(itemPoly));
				n = fread (sim->polymer.items, sizeof(itemPoly), count, chk); if (n != count) error (ECHKREAD);
				break;

			case (BINARY_DESC_FILE):
				// simfiles
				*fpptr = (simfile *) mycalloc (1, sizeof(simfile));
				f = *fpptr;
				n = fread (f, sizeof(simfile), count, chk);	if (n != count) error (ECHKREAD);
				f->next = NULL;
				fpptr = &(f->next);
				break;

			case (BINARY_DESC_HISTOGRAM):

				// histograms
				*hpptr = (simhist *) mycalloc (1, sizeof(simhist));
				h = *hpptr;
				n = fread (h, sizeof(simhist), count, chk);	if (n != count) error (ECHKREAD);

				// read histogram data
				count = h->n;
				HistogramAllocateBins (h);

				n = fread (&desc,	sizeof(int), 1, chk); 		if (n != 1) 	error (ECHKREAD);
				n = fread (&count,   sizeof(int), 1, chk); 		if (n != 1) 	error (ECHKREAD);
				n = fread (h->bin,   sizeof(real), count, chk); 	if (n != count) error (ECHKREAD);

				n = fread (&desc,	sizeof(int), 1, chk); 		if (n != 1) 	error (ECHKREAD);
				n = fread (&count,   sizeof(int), 1, chk); 		if (n != 1) 	error (ECHKREAD);
				n = fread (h->count, sizeof(real), count, chk); 	if (n != count) error (ECHKREAD);

				h->next = NULL;
				hpptr = &(h->next);
				break;

			case (BINARY_DESC_SCENE):

				// scenes
				*spptr = (simscene *) mycalloc (1, sizeof(simscene));
				s = *spptr;
				n = fread (s, sizeof(simscene), count, chk);	if (n != count) error (ECHKREAD);
				s->next = NULL;
				spptr = &(s->next);
				break;

			default:
				error (ECHKREAD);
		}
	}

	// adjust particle pointers
	for (i=0; i<sim->atom.n; i++) {
		if (sim->atom.items[i].next)  sim->atom.items[i].next  = sim->atom.items + (sim->atom.items[i].next  - atomStart);
	}
	for (i=0; i<sim->polymer.n; i++) {
		if (sim->polymer.items[i].p1) sim->polymer.items[i].p1 = sim->atom.items + (sim->polymer.items[i].p1 - atomStart);
	}
}


//================================================================================
void WriteBinary (int desc, int count, int size, void *ptr, FILE *f)
//================================================================================
{
	// Write binary data records to disk. Start with two int fields specifying
	// the type of data and the element count, followed by the data itself.

	int n;

	n = fwrite (&desc, 	sizeof(int), 1, f);  if (n != 1)	 error (ECHKWRITE);
	n = fwrite (&count, sizeof(int), 1, f);  if (n != 1)	 error (ECHKWRITE);
	n = fwrite (ptr, 	size, count, f);	 if (n != count) error (ECHKWRITE);
}


//================================================================================
simfile *GetSimFile (simfile *fptr, char *label)
//================================================================================
{
	// Retrieve the simfile pointer corresponding to a given sim file label.

	int found=0;

	while (!found && fptr) {
		if (strcmp(label, fptr->label) == 0) found = 1;
		else fptr = fptr->next;
	}

	if (!found) {
		error (EFILELIST);
	}

	return (fptr);
}


//================================================================================
FILE *GetSimStream (simfile *fptr, char *label)
//================================================================================
{
	// Retrieve the file stream pointer corresponding to a given sim file label.

	int found=0;

	while (!found && fptr) {
		if (strcmp(label, fptr->label) == 0) found = 1;
		else fptr = fptr->next;
	}

	if (!found) {
		printf ("%s", label);
		error (EFILELIST);
	}

	return (fptr->stream);
}


//================================================================================
char *GetSimFileName (simfile *fptr, char *label)
//================================================================================
{
	// Retrieve the file name corresponding to a given sim file label.

	int found=0;

	while (!found && fptr) {
		if (strcmp(label, fptr->label) == 0) found = 1;
		else fptr = fptr->next;
	}

	if (!found) {
		printf ("%s", label);
		error (EFILELIST);
	}

	return (fptr->name);
}


//================================================================================
simfile *FileNew (simfile *fptr)
//================================================================================
{
	// Adds an output file to the beginning of the linked list of output files
	// pointed to by fptr.

	simfile *new;

	// allocate memory for the new simfile structure
	new = (simfile *) mycalloc (1, sizeof(simfile));

	// insert new element in non-NULL list
	if (fptr) {
		new->next = fptr->next;
		fptr->next = new;
	}

	// return a pointer to the newly created simfile
	return new;
}


//================================================================================
void CloseDataFiles (simptr sim)
//================================================================================
{
	// Closes all the open data files. We could use fcloseall here, but instead
	// we iterate through the file list in order to leave more flexibility in
	// the future if some files need to stay open.

	simfile *fptr;

	// report
	LOG ("Closing all data files\n");
	LOG ("--------------------------------------------------------------\n");

	// close valid streams
	fptr = sim->files;
	while (fptr) {
		if (fptr->stream && fileno(fptr->stream)>0) fclose (fptr->stream);
		fptr = fptr->next;
	}
}
