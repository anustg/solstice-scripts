/* Copyright (C) 2017-2018 |Meso|Star>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>. */

#include "solpp.h"

int
main(int argc, char** argv)
{
  char s[128];
  buf_char_T buf = BUF_NULL;
  FILE* input = stdin;
  FILE* output = NULL;
  char* line = NULL;
  double azim = 0;
  double elev = 0;

  if(argc > 1 && !(input = fopen(argv[1], "r"))) {
    fprintf(stderr, "Could not open the file `%s'.\n", argv[1]);
    return 1;
  }
  while((line = read_line(&buf, input))) {
    if(strncmp(line, "#--- Sun direction:", 19)) {
      CHK(output != NULL);
      fprintf(output, "%s\n", line);
    } else {
      CHK(sscanf(line+19, "%lf %lf (%*f %*f %*f)", &azim, &elev) == 2);
      CHK(snprintf(s, sizeof(s), "%g-%g-paths.vtk", azim, elev) < sizeof(s));
      if(output) fclose(output);
      printf("Writing `%s'\n", s);
      CHK(output = fopen(s, "w"));
    }
  }
  BUF_RELEASE(buf);
  if(output) fclose(output);
  if(input && input!=stdin) fclose(input);
  return 0;
}
