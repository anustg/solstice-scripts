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
  char* line = NULL;
  double azim = -1;
  double elev = -1;

  if(argc > 1 && !(input = fopen(argv[1], "r"))) {
    fprintf(stderr, "Could not open the file `%s'.\n", argv[1]);
    return 1;
  }

  line = read_line(&buf, input);
  while(line) {
    if(!strncmp(line, "#--- Sun direction:", 19)) {
      /* Get the solar direction */
      CHK(sscanf(line+19, "%lf %lf (%*f %*f %*f)", &azim, &elev)==2);
      line = read_line(&buf, input);
    } else if(!strncmp(line, "# vtk", 5)) {
      char* header = NULL;
      char* rcv_name = NULL;
      FILE* output;

      CHK(azim >= 0 && elev >= 0);

      CHK(header = strdup(line)); /* Duplicate the current line */
      CHK(line = read_line(&buf, input));
      CHK(rcv_name = strdup(line)); /* Duplicate the line of the receiver name */

      /* Create the name of the destination file */
      CHK(snprintf(s, sizeof(s), "%g-%g-%s.vtk", azim, elev, rcv_name) < sizeof(s));
      printf("Writing `%s'\n", s);
      CHK(output = fopen(s, "w"));

      /* Write the map data into `output' */
      fprintf(output, "%s\n", header);
      fprintf(output, "%s\n", rcv_name);
      while((line = read_line(&buf, input)) && line[0] != '#') {
        fprintf(output, "%s\n", line);
      }

      /* Clean up temporary variable and close the destination file */
      fclose(output);
      free(header);
      free(rcv_name);
    } else {
      line = read_line(&buf, input);
    }
  }

  BUF_RELEASE(buf);
  if(input && input!=stdin) fclose(input);
  return 0;
}

