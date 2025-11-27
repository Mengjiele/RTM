#include "par.cuh"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdbool.h>
#define MAX_PARAMS 100
#define MAX_LINE_LENGTH 256

typedef struct {
    char key[64];
    char value[128];
} Parameter;

static Parameter params[MAX_PARAMS];
static int param_count = 0;

static void load_parameters(const char *filename)
{
    FILE *file = fopen(filename, "r");
    if (!file) err("Cannot open parameter file");

    char line[MAX_LINE_LENGTH];
    while (fgets(line, MAX_LINE_LENGTH, file))
	{
		 if (param_count >= MAX_PARAMS)
		 {
            fprintf(stderr, "Warning: Too many parameters (max=%d)\n", MAX_PARAMS);
            break;
        }
		
        if (line[0] == '#' || line[0] == '\n') continue;
        
        char *equal_pos = strchr(line, '=');
        if (!equal_pos) continue;
        
        *equal_pos = '\0';
        char *key = line;
        char *value = equal_pos + 1;
        
        value[strcspn(value, "\n")] = '\0';
        
        if (param_count < MAX_PARAMS)
		{
            strncpy(params[param_count].key, key, 63);
            strncpy(params[param_count].value, value, 127);
            param_count++;
        }
    }
    fclose(file);
}

void parse_command_line(int argc, char **argv)
{
    for (int i = 1; i < argc; i++) {
        if (strncmp(argv[i], "par=", 4) == 0)
		{
            load_parameters(argv[i] + 4);
        }
    }
    if (param_count == 0) err("No parameter file specified (use par=filename.par)");
}

int getparint(const char *key, int *val)
{
    for (int i = 0; i < param_count; i++)
	{
        if (strcmp(params[i].key, key) == 0)
		{
            *val = atoi(params[i].value);
            return 1;
        }
    }
    return 0;
}


int getparfloat(const char *key, float *val)
{
    for (int i = 0; i < param_count; i++)
	{
        if (strcmp(params[i].key, key) == 0)
		{
            *val = atof(params[i].value);
            return 1;
        }
    }
    return 0;
}


int getparstring(const char *key, char *val)
{
    for(int i = 0; i < param_count; i++)
	{
        if(strcmp(params[i].key, key) == 0)
		{
            strncpy(val, params[i].value, 127);
            val[127] = '\0';
            return 1;
        }
    }
    return 0;
}



int getparbool(const char *key, bool *val)
{
    char strval[16] = {0};
    if (!getparstring(key, strval))
	{
        return 0;
    }

    *val = (strncasecmp(strval, "yes", 3) == 0 || 
           strncasecmp(strval, "true", 4) == 0 ||
           atoi(strval) != 0);
    return 1;
}

void err(const char *msg)
{
    fprintf(stderr, "ERROR: %s\n", msg);
    exit(EXIT_FAILURE);
}