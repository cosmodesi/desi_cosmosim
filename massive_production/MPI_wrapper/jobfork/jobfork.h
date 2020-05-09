/*******************************************************************************
* 
* jobfork: C tool for running multiple jobs in parallel.
*
* Github repository:
*       https://github.com/cheng-zhao/jobfork
*
* Copyright (c) 2020 Cheng Zhao <zhaocheng03@gmail.com>
* 
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
* 
* The above copyright notice and this permission notice shall be included in all
* copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
* SOFTWARE.
*
*******************************************************************************/

#ifndef __JOBFORK_H__
#define __JOBFORK_H__

#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <signal.h>
#include <sys/types.h>

#ifdef CMD_MPI
#define FLUSH_STDOUT
#undef CMD_OMP
#endif

#ifndef _PATH_BSHELL
#define _PATH_BSHELL "/bin/sh"
#endif

#define CMD_BUF 2048            /* maximum number of characters for commands */
#define COMMENT '#'             /* comment symbol for the job list file */
#define JOB_START 1
#define JOB_FAIL  2
#define JOB_DONE  3

volatile sig_atomic_t term;     /* flag for signal termination */

/******************************************************************************
  Definition of error codes.
******************************************************************************/
#define ERR_MEMORY      101     /* failed to allocate memory      */
#define ERR_FILE        102     /* failed to read file            */
#define ERR_PIPE        103     /* failed to allocate pipe        */
#define ERR_FORK        104     /* failed to create a process     */
#define ERR_REDIR       105     /* failed to redirect I/O         */
#define ERR_EXEC        106     /* failed to execute the job      */
#define ERR_STRING      107     /* string length exceeding limits */
#define ERR_ARG         108     /* invalid argument               */
#define ERR_CMD         109     /* invalid command                */
#define ERR_OTHER       199     /* unknown errors                 */

/******************************************************************************
  Definitions for printing messages.
******************************************************************************/
#define MSG_ERR(...)            \
  fprintf(stderr, "\x1B[31;1mError:\x1B[0m " __VA_ARGS__)
#define MSG_CHILD_STDOUT(...)   \
  printf("\x1B[32;1m[%d-%d]\x1B[0m %s", __VA_ARGS__)
#define MSG_CHILD_STDERR(...)   \
  fprintf(stderr, "\x1B[31;1m<%d-%d>\x1B[0m %s", __VA_ARGS__)

/******************************************************************************
  Definition of data types.
******************************************************************************/
typedef struct {
  pid_t pid;
  FILE *out;
  FILE *err;
} CHILD_INFO;

struct cmd_status {
  char fname_rst[CMD_BUF];              /* restart file for failed jobs */
  int num;
  int len;
  char *cmd;
  char *status;
} cstat;

/******************************************************************************
  Definition of functions.
******************************************************************************/
int read_jobs(const char *);

void save_jobs(void);

int create_child(const char *, CHILD_INFO *);

int close_child(CHILD_INFO *);

void terminate(int);

void mpi_manager(const int, int *);

void mpi_worker(char *, CHILD_INFO *);

#endif
