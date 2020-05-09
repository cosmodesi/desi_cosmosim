#include "jobfork.h"
#include <ctype.h>

#ifdef CMD_MPI
#include <mpi.h>
#else
#ifdef CMD_OMP
#include <omp.h>
#endif
#endif

int main(int argc, char *argv[]) {
  char *cmd;
  int ncmd, ret;
  CHILD_INFO ci;

  if (argc != 2) {
    fprintf(stderr, "Usage: %s job_list\n"
        "joblist: a text file with each line being a job (command)\n",
        argv[0]);
    return ERR_ARG;
  }
  term = 0;

#ifdef CMD_MPI
  int myrank, tasknum;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &tasknum);

  if (tasknum < 2) {
    MSG_ERR("unable to distribute jobs with %d MPI tasks.\n", tasknum);
    MPI_Finalize();
    return ERR_OTHER;
  }

  if (myrank == 0) {            /* manager */
#endif

  if ((ret = read_jobs(argv[1]))) return ret;
  ncmd = cstat.num;
  printf("%d jobs are found in the list file.\n", ncmd);

  ret = snprintf(cstat.fname_rst, CMD_BUF, "%s.rst", argv[1]);
  if (ret < 0 || ret >= CMD_BUF) {
    MSG_ERR("the filename is too long to create the restart file.\n");
    return ERR_STRING;
  }

  cstat.status = calloc(ncmd, sizeof(char));
  if (!(cstat.status)) {
    MSG_ERR("failed to allocate memory for recording the job status.\n");
    return ERR_MEMORY;
  }

#ifdef CMD_MPI
  printf("Parallelising jobs with MPI: %d tasks, %d workers.\n",
      tasknum, tasknum - 1);
  }
  MPI_Bcast(&cstat.len, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&cstat.num, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

#ifdef CMD_MPI
  if (myrank == 0) {            /* signal handling only by manager */
#endif
    if (atexit(save_jobs)) {
      MSG_ERR("unable to register exit functions.\n"
          "Restart file may not be created.\n");
    }

    struct sigaction sa;
    sa.sa_handler = &terminate;
    sa.sa_flags = SA_RESTART;
    sigfillset(&sa.sa_mask);
    if (sigaction(SIGHUP, &sa, NULL) == -1 || 
        sigaction(SIGINT, &sa, NULL) == -1 ||
        sigaction(SIGTERM, &sa, NULL) == -1 ||
        sigaction(SIGQUIT, &sa, NULL) == -1) {
      MSG_ERR("unable to catch signals.\n"
          "Restart file may not be created.\n");
    }

#ifdef CMD_MPI
    /* jobs distributed by manager */
    int *nsent;
    nsent = calloc(tasknum - 1, sizeof(int));
    if (!nsent) {
      MSG_ERR("failed to allocate memory for the manager.\n");
      MPI_Abort(MPI_COMM_WORLD, ERR_MEMORY);
      return ERR_MEMORY;
    }
    mpi_manager(tasknum - 1, nsent);
    free(nsent);
  }
  else {                        /* workers */
    cmd = calloc(cstat.len, sizeof(char));
    if (!cmd) {
      MSG_ERR("failed to allocate memory for the workers.\n");
      MPI_Abort(MPI_COMM_WORLD, ERR_MEMORY);
      return ERR_MEMORY;
    }
    mpi_worker(cmd, &ci);
    free(cmd);
  }

  MPI_Finalize();

#else

#ifdef CMD_OMP

  printf("Parallelising jobs with OpenMP: %d threads.\n",
      omp_get_max_threads());
  int id = 0;
  char line[CMD_BUF];

#pragma omp parallel for private(ci,line,cmd) firstprivate(id) schedule(dynamic)
  for (int i = 0; i < ncmd; i++) {
    cmd = cstat.cmd + i * cstat.len;
    printf("-> Allocating command to thread %d (job index: %d):\n   %s\n",
        omp_get_thread_num(), id, cmd);
    cstat.status[i] = JOB_START;
    if ((ret = create_child(cmd, &ci))) {
      MSG_ERR("failed to execute command `%s'.\n", cmd);
      cstat.status[i] = JOB_FAIL;
    }
    else {
      memset(line, 0, CMD_BUF);
      while (fgets(line, CMD_BUF, ci.out) != NULL) {
        MSG_CHILD_STDOUT(omp_get_thread_num(), id, line);
        memset(line, 0, CMD_BUF);
      }

      while (fgets(line, CMD_BUF, ci.err) != NULL) {
        MSG_CHILD_STDERR(omp_get_thread_num(), id, line);
        memset(line, 0, CMD_BUF);
      }

      if ((ret = close_child(&ci))) {
        MSG_ERR("unable to finish command `%s'.\n", cmd);
        cstat.status[i] = JOB_FAIL;
      }
      else cstat.status[i] = JOB_DONE;
    }

    id++;
  }

#endif

#endif

  return 0;
}


void terminate(int sig) {
  save_jobs();
}

