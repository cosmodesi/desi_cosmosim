#include "jobfork.h"
#include <sys/wait.h>

#define PIPE_READ 0
#define PIPE_WRITE 1

extern char **environ;

int create_child(const char *cmd, CHILD_INFO *ci) {
  int pipe_stdout[2];
  int pipe_stderr[2];
  int ret;
  pid_t pid;
  char *argp[] = {"sh", "-c", NULL, NULL};

  if (pipe(pipe_stdout) < 0) {
    close(pipe_stdout[PIPE_READ]);
    close(pipe_stdout[PIPE_WRITE]);
    MSG_ERR("failed to allocate pipe for child stdout on job:\n"
        "    %s\n", cmd);
    return ERR_PIPE;
  }
  if (pipe(pipe_stderr) < 0) {
    close(pipe_stderr[PIPE_READ]);
    close(pipe_stderr[PIPE_WRITE]);
    MSG_ERR("failed to allocate pipe for child stderr on job:\n"
        "    %s\n", cmd);
    return ERR_PIPE;
  }

  pid = fork();
  if (pid < 0) {
    close(pipe_stdout[PIPE_READ]);
    close(pipe_stdout[PIPE_WRITE]);
    close(pipe_stderr[PIPE_READ]);
    close(pipe_stderr[PIPE_WRITE]);
    MSG_ERR("failed to create a new process for job:\n"
        "    %s\n", cmd);
    return ERR_FORK;
  }
  else if (pid == 0) {          /* child process */
    if (dup2(pipe_stdout[PIPE_WRITE], STDOUT_FILENO) == -1) {
      MSG_ERR("failed to redirect stdout for job:\n"
          "    %s\n", cmd);
      return ERR_REDIR;
    }
    if (dup2(pipe_stderr[PIPE_WRITE], STDERR_FILENO) == -1) {
      MSG_ERR("failed to redirect stderr for job:\n"
          "    %s\n", cmd);
      return ERR_REDIR;
    }
    close(pipe_stdout[PIPE_READ]);
    close(pipe_stdout[PIPE_WRITE]);
    close(pipe_stderr[PIPE_READ]);
    close(pipe_stderr[PIPE_WRITE]);

    /* run child process */
    argp[2] = (char *) cmd;
    ret = execve(_PATH_BSHELL, argp, environ);
    exit(ret);
  }

  /* pid > 0, parent process */
  close(pipe_stdout[PIPE_WRITE]);
  close(pipe_stderr[PIPE_WRITE]);

  ci->out = fdopen(pipe_stdout[PIPE_READ], "r");
  ci->err = fdopen(pipe_stderr[PIPE_READ], "r");
  ci->pid = pid;

  return 0;
}


int close_child(CHILD_INFO *ci) {
  int pstat = 1;
  pid_t pid;

  if (ci == NULL) return ERR_OTHER;
  fclose(ci->out);
  fclose(ci->err);
  do {
    pid = waitpid(ci->pid, &pstat, 0);
  } while (pid == -1);

  return (pid == -1 ? -1 : pstat);
}

