#include "jobfork.h"
#include <ctype.h>
#include <fcntl.h>

#ifdef CMD_MPI
#include <mpi.h>
#endif

int read_jobs(const char *fname) {
  FILE *fp;
  char line[CMD_BUF];
  char *buf;
  int i, n, maxlen;

  if (!(fp = fopen(fname, "r"))) {
    MSG_ERR("cannot open the job list file `%s'.\n", fname);
    return ERR_FILE;
  }

  /* count the number and maximum length of commands */
  n = maxlen = 0;
  memset(line, 0, CMD_BUF);
  while (fgets(line, CMD_BUF, fp) != NULL) {
    sscanf(line, "%*[ |\t]%[^\n]", line);       /* remove leading spaces */
    if (line[0] != COMMENT && isgraph(line[0])) {
      buf = memchr(line, '\n', CMD_BUF);
      if (buf == NULL) {
        MSG_ERR("command length exceeds CMD_BUF.\n");
        return ERR_STRING;
      }

      i = buf - line + 1;
      if (maxlen < i) maxlen = i;
      n++;
    }
    memset(line, 0, CMD_BUF);
  }

  if (n < 1) {
    MSG_ERR("job not found.\n");
    return ERR_FILE;
  }

  cstat.num = n;
  cstat.len = maxlen;
  cstat.cmd = calloc((size_t) n * maxlen, sizeof(char));
  if (!(cstat.cmd)) {
    MSG_ERR("failed to allocate memory for the commands.\n");
    return ERR_MEMORY;
  }

  /* read commands into array */
  fseek(fp, 0, SEEK_SET);
  n = 0;
  while (fgets(line, CMD_BUF, fp) != NULL) {
    sscanf(line, "%*[ |\t]%[^\n]", line);       /* remove leading spaces */
    if (line[0] != COMMENT && isgraph(line[0])) {
      line[strcspn(line, "\r\n")] = '\0';
      buf = cstat.cmd + n * maxlen;
      strncpy(buf, line, maxlen);
      n++;
    }
    memset(line, 0, CMD_BUF);
  }

  fclose(fp);
  return 0;
}


void save_jobs(void) {
  if (term == 1)  return;       /* this function should only run for once */
  else term = 1;

  int i, cnt, logfile;
  char *cmd;

  for (cnt = 0, i = 0; i < cstat.num; i++) {
    if (cstat.status[i] != JOB_DONE) cnt++;
  }
  if (cnt == 0) return;         /* no failed job */

  if (!(cstat.status && cstat.cmd)) return;

  logfile = open(cstat.fname_rst, O_WRONLY | O_CREAT | O_TRUNC, 0644);
  write(logfile, "# Unfinished jobs\n", 18);

  for (i = 0; i < cstat.num; i++) {
    if (cstat.status[i] != JOB_DONE) {
      cmd = cstat.cmd + i * cstat.len;
      write(logfile, cmd, strlen(cmd));
      write(logfile, "\n", 1);
    }
  }

  write(STDOUT_FILENO, "Restart file created:\n  ", 24);
  write(STDOUT_FILENO, cstat.fname_rst, strlen(cstat.fname_rst));
  write(STDOUT_FILENO, "\n", 1);
  close(logfile);
}

