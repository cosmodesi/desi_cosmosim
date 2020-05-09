#include "jobfork.h"
#include <mpi.h>

void mpi_manager(const int nworker, int *nsent) {
  int i, n, numsent, numrecv, sender, irecv, job_status;
  char *cmd;
  MPI_Status status;

  /* send one job to each worker */
  n = (nworker < cstat.num) ? nworker : cstat.num;
  numsent = numrecv = 0;
  for (i = 0; i < n; i++) {
    cmd = cstat.cmd + numsent * cstat.len;
    printf("-> Allocating command to task %d (job index: %d):\n   %s\n",
        i + 1, nsent[i], cmd);
#ifdef FLUSH_STDOUT
    fflush(stdout);
#endif
    MPI_Send(nsent + i, 1, MPI_INT, i + 1, 0, MPI_COMM_WORLD);
    MPI_Send(cmd, cstat.len, MPI_CHAR, i + 1, numsent, MPI_COMM_WORLD);
    numsent++;
    nsent[i] += 1;
  }

  /* send jobs to free workers */
  for (i = numsent; i < cstat.num; i++) {
    MPI_Recv(&job_status, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
        MPI_COMM_WORLD, &status);

    sender = status.MPI_SOURCE;
    irecv = status.MPI_TAG;
    cstat.status[irecv] = job_status;
    numrecv++;

    cmd = cstat.cmd + i * cstat.len;
    printf("-> Allocating command to task %d (job index: %d):\n   %s\n",
        sender, nsent[sender - 1], cmd);
#ifdef FLUSH_STDOUT
    fflush(stdout);
#endif
    MPI_Send(nsent + sender - 1, 1, MPI_INT, sender, 0, MPI_COMM_WORLD);
    MPI_Send(cmd, cstat.len, MPI_CHAR, sender, numsent, MPI_COMM_WORLD);
    numsent++;
    nsent[sender - 1] += 1;
  }

  /* no more jobs to be sent */
  i = -1;
  while (numrecv < cstat.num) {
    MPI_Recv(&job_status, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG,
        MPI_COMM_WORLD, &status);

    sender = status.MPI_SOURCE;
    irecv = status.MPI_TAG;
    cstat.status[irecv] = job_status;

    numrecv++;
    MPI_Send(&i, 1, MPI_INT, sender, 0, MPI_COMM_WORLD);
  }
}


void mpi_worker(char *cmd, CHILD_INFO *ci) {
  int idx = 0, cnt, job_status, ret, myrank;
  char line[CMD_BUF];
  MPI_Status status;

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Recv(&cnt, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

  while (cnt >= 0) {
    MPI_Recv(cmd, cstat.len, MPI_CHAR, 0, MPI_ANY_TAG,
        MPI_COMM_WORLD, &status);
    idx = status.MPI_TAG;

    job_status = JOB_START;
    if ((ret = create_child(cmd, ci))) {
      MSG_ERR("failed to execute command `%s'.\n", cmd);
      job_status = JOB_FAIL;
    }
    else {
      memset(line, 0, CMD_BUF);
      while (fgets(line, CMD_BUF, ci->out) != NULL) {
        MSG_CHILD_STDOUT(myrank, cnt, line);
#ifdef FLUSH_STDOUT
        fflush(stdout);
#endif
        memset(line, 0, CMD_BUF);
      }

      while (fgets(line, CMD_BUF, ci->err) != NULL) {
        MSG_CHILD_STDERR(myrank, cnt, line);
#ifdef FLUSH_STDOUT
        fflush(stderr);
#endif
        memset(line, 0, CMD_BUF);
      }

      if ((ret = close_child(ci))) {
        MSG_ERR("unable to finish command `%s'.\n", cmd);
#ifdef FLUSH_STDOUT
        fflush(stdout);
#endif
        job_status = JOB_FAIL;
      }
      else job_status = JOB_DONE;
    }

    MPI_Send(&job_status, 1, MPI_INT, 0, idx, MPI_COMM_WORLD);
    MPI_Recv(&cnt, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
  }
}

