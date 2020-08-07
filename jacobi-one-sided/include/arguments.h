#ifndef ARGUMENTS_H
#define ARGUMENTS_H

#include "data.h"

#define ARGUMENT_NAMES_MAX_LENGTH 64
#define ARGUMENT_DESCRIPTION_MAX_LENGTH 256
#define NO_OPTION ""
#define ARGUMENTS_SIZE 5
#define ERROR_UNKNOWN_ARGUMENT -1
#define ERROR_INVALID_OPTION -2
#define ERROR_PARSING_ARGUMENTS -1

typedef struct
{
	int index;
	char full_name[ARGUMENT_NAMES_MAX_LENGTH];
	char short_name[ARGUMENT_NAMES_MAX_LENGTH];
	char option_name[ARGUMENT_NAMES_MAX_LENGTH];
	char description[ARGUMENT_DESCRIPTION_MAX_LENGTH];
	char* option_value;
} argument_t;

int parse_command_line_arguments(int argc, char* argv[], instance_t* instance);
void parse_argument(char* arg);
void print_help();

#endif // !ARGUMENTS_H