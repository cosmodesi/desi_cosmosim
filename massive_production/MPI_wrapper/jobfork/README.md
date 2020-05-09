# jobfork

![GitHub](https://img.shields.io/github/license/cheng-zhao/libcfg.svg)

## Table of Contents

-   [Introduction](#introduction)
-   [Compilation](#compilation)
-   [Getting started with examples](#getting-started-with-examples)
    -   [Job list file](#job-list-file)
    -   [OpenMP schedulerr](#openmp-scheduler)
    -   [MPI scheduler](#mpi-scheduler)
    -   [Standard output/error](#standard-outputerror)
    -   [Restart file](#restart-file)

## Introduction

jobfork is a simple tool written in C, for running a list of commands in parallel, with either OpenMP or MPI schedulers. It is designed to run multiple jobs simultaneously (but independently) with the granted resources, and launch a pending task in the list as soon as any job has finished.

Jobs that are already parallelised with OpenMP or MPI cannot be run by jobfork compiled the same parallelisation scheme. However, as a common practice, OpenMP-parallelised jobs can be called by jobfork with the MPI scheduler.

jobfork captures all the standard outputs and errors of the jobs, and prepend identifiers with different symbols and ANSI colours. Jobs that return non-zero states and unfinished jobs are saved to a restart file on the termination of jobfork.

Apart from the reliance on `omp.h` and `mpi.h`, which are required for the OpenMP and MPI schedulers respectively, jobfork is compliant with the ISO C99 and IEEE POSIX.1-2008 standards. It is written by Cheng Zhao (&#36213;&#25104;), and is distributed under the MIT license (see [LICENSE.txt](LICENSE.txt) for the details).

If you use this tool in research that results in publications, please consider citing the following paper:
```
Zhao et al. in preparation.
```

<small>[\[TOC\]](#table-of-contents)</small>

## Compilation

Depending on the scheduler, jobfork has different dependencies in addition to the C99 standard library and the C POSIX library, i.e., `omp.h` for OpenMP and `mpi.h` for MPI, as well as the C compiler that supports the corresponding parallelisation libraries. In particular, for the MPI scheduler, a MPI library (such as [Open MPI](https://www.open-mpi.org/)) has to be configured for some supercomputers (for instance, with the `module load` command).

By default, jobfork can be compiled with both schedulers using the command
```bash
make
```
Upon successful completion of the compilation process, two executables `jobfork_omp` and `jobfork_mpi` are generated, for the OpenMP and MPI schedulers, respectively.

Alternatively, the two components can be compiled individually with the following commands
```bash
make omp    # for the OpenMP scheduler
make mpi    # for the MPI scheduler
```
In this case only the corresponding executable is generated for each of the compilation command.

<small>[\[TOC\]](#table-of-contents)</small>

## Getting started with examples

### Job list file

Before calling the jobfork executables, one has to prepare a job list file, which is a plain text file containing all the jobs to be run. Each line of the file should be one valid shell command. Lines starting with `#` (which can be reset in [jobfork.h](jobfork.h#L55)) are treated as comments, and ignored by the scheduler.

Note that currently the job list file must be a real file on disk, and should be passed by pipe. In addition, the number of characters for each command should not exceed 2048 (defined in [jobfork.h](jobfork.h#L54)).

As an example, we create a job list file called `joblist.txt`, with the contents as follows:
```console
$ cat joblist.txt
echo "Hello World!"
uname
non-existent-command
```

There are three commands in total, in which the last one is not valid and should fail eventually.

<small>[\[TOC\]](#table-of-contents)</small>

### OpenMP scheduler

To run jobfork with OpenMP, the number of threads to be used should be specified via the environment variable `OMP_NUM_THREADS`. Then the `jobfork_omp` executable can be simply called followed by the job list file.

Below is an example for running the jobs in `joblist.txt` in parallel with two OpenMP threads:
```console
$ export OMP_NUM_THREADS=2
$ jobfork_omp joblist.txt
3 jobs are found in the list file.
Parallelising jobs with OpenMP: 2 threads.
-> Allocating command to thread 0 (job index: 0):
   echo "Hello World!"
-> Allocating command to thread 1 (job index: 0):
   uname
[0-0] Hello World!
-> Allocating command to thread 0 (job index: 1):
   non-existent-command
[1-0] Linux
<0-1> sh: non-existent-command: command not found
Error: unable to finish command `non-existent-command'.
Restart file created:
  joblist.txt.rst
```

Here, the first two commands are assigned to the two OpenMP threads with thread IDs being `0` and `1`, respectively. They are the first commands run on their threads, hence are both identified by the job index `0`. Then the thread `0` finishes its job first, and is allocated the third command `non-existent-command`, with the job index incremented by 1.

<small>[\[TOC\]](#table-of-contents)</small>

### MPI scheduler

The `jobfork_mpi` executable can be called in the standard way for running MPI programs, with a MPI process manager such as `mpirun`. The number of MPI tasks can then be specified via command line options of the process manager.
However, not all of the MPI tasks are used to run the jobs in parallel. There is one task acting as the "manager", which maintains the job pool and distribute jobs to the other tasks that are dubbed "workers". Thus, the actual number of jobs running in parallel is always one less than the total number of MPI tasks.

An example for running the jobs in `joblist.txt` using two MPI tasks is shown below:
```console
$ mpirun -n 2 ./jobfork_mpi joblist.txt
3 jobs are found in the list file.
Parallelising jobs with MPI: 2 tasks, 1 workers.
-> Allocating command to task 1 (job index: 0):
   echo "Hello World!"
[1-0] Hello World!
-> Allocating command to task 1 (job index: 1):
   uname
[1-1] Linux
-> Allocating command to task 1 (job index: 2):
   non-existent-command
<1-2> sh: non-existent-command: command not found
Error: unable to finish command `non-existent-command'.
Restart file created:
  joblist.txt.rst
```

It can be seen that with only two MPI tasks, the jobs are essentially run in sequence, as there is only one worker.

<small>[\[TOC\]](#table-of-contents)</small>

### Standard output/error

The standard outputs and errors of each jobs is reprocessed by jobfork, and prepended a unique identifier, which consists of two numbers connected by the `-` symbol. In particular, the two numbers are the corresponding thread ID and job index, respectively. Furthermore, the identifiers for standard outputs are surrounded by pairs of square brackets (`[]`) with green ANSI colour, while for standard errors they are red angle brackets (`<>`).

<small>[\[TOC\]](#table-of-contents)</small>

### Restart file

If not all of the jobs finish successfully (recognised by 0 exit codes) on the termination by jobfork, even if jobfork is aborted abnormally (such as being cancelled by the user with the `CTRL`-`C` keys), a restart file is created with all the failed or unfinished commands. The restart file is in the same format as the job list file, hence can be run with jobfork later.

The name of the restart file is hard-coded, with the suffix `.rst` appended to the name of the job list file. If the file already exists, then jobfork will overwrite it silently.

In our example job list file, there is a nonexistent command `non-existent-command`, which cannot be executed by shell. As a failed job, it is recorded in the restart file for both the OpenMP and MPI schedulers:
```console
$ cat joblist.txt.rst
# Unfinished jobs
non-existent-command
```

<small>[\[TOC\]](#table-of-contents)</small>
