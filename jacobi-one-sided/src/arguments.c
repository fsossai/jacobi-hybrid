#include <string.h>
#include <stdlib.h>

#include "arguments.h"

static argument_t arguments[ARGUMENTS_SIZE] =
{
	{ 0,"no-shared-memory",			"m",	NO_OPTION,		"Choose whether to use the shared memory."	},
	{ 1,"split-direction",			"s",	"DIRECTION",	"Direction along which to split the workload "
															"in the shared memory. Available values are "
															"X, Y or Z. Default value is X." },
	{ 2,"cart-no-split-direction",	"n",	"DIRECTION",	"Select along which direction the cartesian "
															"topology should be flattened. Available "
															"values are X, Y or Z." },
	{ 3,"heads-per-shared-region",	"H",	"NUM",			"Set how many 'head' processes will be "
															"instantiated per each physical shared "
															"region of memory. Default value is 1." },
	{ 4,"input-instance",			"i",	"FILENAME",		"Input file with domain sizes and parameters "
															"one per line: size x, size y, size z, alpha, "
															"relaxation, tolerance, max iterations." },
	{ 5,"help",						"h",	NO_OPTION,		"Show this help message." }
};

int parse_command_line_arguments(int argc, char* argv[], instance_t* instance)
{
	instance->use_shared_memory = DEFAULT_USE_SHARED_MEMORY;
	instance->local_subdomain_split_direction = DEFAULT_SHARED_SPLIT_DIRECTION;
	instance->heads_per_shared_region = DEFAULT_HEADS_PER_SHARED_REGION;
	memset(instance->cart_splits, 0x00, DOMAIN_DIM * sizeof(int));

	int n;
	for (int i = 1; i < argc; i++)
	{
		argument_t argument;
		parse_argument(argv[i], &argument);

		switch (argument.index)
		{
		case 0:
			instance->use_shared_memory = 0;
			break;
		case 1:
			if (!strcmp(argument.option_value, "X") || !strcmp(argument.option_value, "x"))
				instance->local_subdomain_split_direction = 0;
			else if (!strcmp(argument.option_value, "Y") || !strcmp(argument.option_value, "y"))
				instance->local_subdomain_split_direction = 1;
			else if (!strcmp(argument.option_value, "Z") || !strcmp(argument.option_value, "z"))
				instance->local_subdomain_split_direction = 2;
			else
			{
				fprintf(stderr, "ERROR: Invalid or missing option for argument '%s'\n", argument.full_name);
				return ERROR_PARSING_ARGUMENTS;
			}
			break;
		case 2:
			if (!strcmp(argument.option_value, "X") || !strcmp(argument.option_value, "x"))
				instance->cart_splits[0] = 1;
			else if (!strcmp(argument.option_value, "Y") || !strcmp(argument.option_value, "y"))
				instance->cart_splits[1] = 1;
			else if (!strcmp(argument.option_value, "Z") || !strcmp(argument.option_value, "z"))
				instance->cart_splits[2] = 1;
			else
			{
				fprintf(stderr, "ERROR: Invalid or missing option for argument '%s'\n", argument.full_name);
				return ERROR_PARSING_ARGUMENTS;
			}
			break;
		case 3:
			n = atoi(argument.option_value);
			if (n > 0)
				instance->heads_per_shared_region = n;
			else
			{
				fprintf(stderr, "ERROR: Invalid or missing option for argument '%s'\n", argument.full_name);
				return ERROR_PARSING_ARGUMENTS;
			}
			break;
		case 4:
			strcpy(instance->input_file_name, argument.option_value);
			instance->input_stream = fopen(instance->input_file_name, "r");
			if (instance->input_stream == NULL)
			{
				fprintf(stderr, "ERROR: Cannot open '%s'\n", instance->input_file_name);
				return ERROR_PARSING_ARGUMENTS;
			}
			break;
		case 5:
			print_help();
			return HELP_MESSAGE;
		case ERROR_UNKNOWN_ARGUMENT:
			fprintf(stderr, "ERROR: Unknown argument '%s'\n", argument.full_name);
			print_help();
			return ERROR_UNKNOWN_ARGUMENT;
		case ERROR_INVALID_OPTION:
			fprintf(stderr, "ERROR: Invalid or missing option for argument '%s'\n", argument.full_name);
			return ERROR_PARSING_ARGUMENTS;
		default:
			return ERROR_PARSING_ARGUMENTS;
		}
	}
	return 0;
}

void parse_argument(char* arg_str, argument_t* argument)
{
	size_t arg_str_length = strlen(arg_str);
	char* token = strtok(arg_str, "=");
	strcpy(argument->full_name, token);
	if (token == NULL || strlen(token) < 2)
	{
		argument->index = ERROR_UNKNOWN_ARGUMENT;
		return;
	}

	for (int i = 0; i < ARGUMENTS_SIZE; i++)
	{
		if (!strcmp(&token[1], arguments[i].short_name) ||
			!strcmp(&token[2], arguments[i].full_name))
		{
			argument->index = i;
			if (!strcmp(arguments[i].option_name, NO_OPTION))
				return;

			char* option = &arg_str[strlen(token) + 1];
			if (arg_str_length - strlen(token) < 2 || strlen(option) == 0)
			{
				argument->index = ERROR_INVALID_OPTION;
				return;
			}

			argument->option_value = option;
			return;
		}
	}
	argument->index = ERROR_UNKNOWN_ARGUMENT;
}

void print_help()
{
	printf("Solving the Helmoltz equation on a %iD grid.\n", DOMAIN_DIM);
	printf("\nUsage: jacobi ");
	for (int i = 0; i < ARGUMENTS_SIZE; i++)
	{
		printf("[-%s%s%s] ", arguments[i].short_name,
			(!strcmp(arguments[i].option_name, NO_OPTION)) ? "" : "=",
			(!strcmp(arguments[i].option_name, NO_OPTION)) ? "" : arguments[i].option_name);
	}

	printf("\n\nIf no input file is specified, the instance data will be requested in the standard input (stdin).\n");
	printf("List of options:\n");

	for (int i = 0; i < ARGUMENTS_SIZE; i++)
	{
		printf("\n -%s, --%s%s%s\n\t%s\n",
			arguments[i].short_name,
			arguments[i].full_name,
			(!strcmp(arguments[i].option_name, NO_OPTION)) ? "" : "=",
			(!strcmp(arguments[i].option_name, NO_OPTION)) ? "" : arguments[i].option_name,
			arguments[i].description);
	}
	printf("\n");
	printf("Author\t\t: Federico Sossai, federico.sossai@gmail.com\n");
	printf("Source code\t: https://github.com/fsossai/jacobi-hybrid\n");
}