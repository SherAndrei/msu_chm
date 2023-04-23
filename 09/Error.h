#pragma once

#include <stdio.h>

enum Error { Success, IncorrectUsage, InputError, LogicError, NotEnoughMemory };

inline int ExplainError(int error, const char* hint) {
  if (error == Success)
    return error;

  switch (error)
  {
  case IncorrectUsage:
    fprintf(stderr, "error: incorrect usage\n");
    break;
  case InputError:
    fprintf(stderr, "error: parsing %s from input\n", hint);
    break;
  case LogicError:
    fprintf(stderr, "error: %s", hint);
    break;
  case NotEnoughMemory:
    fprintf(stderr, "error: not enough memory for %s\n", hint);
    break;
  }
  return error;
}
