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

static inline void
print_rcv_side(FILE* output, const enum side side, struct rcv* rcv)
{
  #define W(Name) rcv->Name[side].E, rcv->Name[side].SE
  fprintf(output, "                 |              [Incoming]                ");
  fprintf(output, "              [Absorbed]               \n");
  fprintf(output, "            Flux | %16g +/- %-16g | %16g +/- %-16g\n",
    W(in.flux), W(abs.flux));
  fprintf(output, "   Material loss | %16g +/- %-16g | %16g +/- %-16g\n",
    W(in.flux_mat_loss), W(abs.flux_mat_loss));
  fprintf(output, "Atmospheric loss | %16g +/- %-16g | %16g +/- %-16g\n",
    W(in.flux_atm_loss), W(abs.flux_atm_loss));
  fprintf(output, "No Material loss | %16g +/- %-16g | %16g +/- %-16g\n",
    W(in.flux_no_mat_loss), W(abs.flux_no_mat_loss));
  fprintf(output, "  No Atmos. loss | %16g +/- %-16g | %16g +/- %-16g\n",
    W(in.flux_no_atm_loss), W(abs.flux_no_atm_loss));
  fprintf(output, "                 |\n");
  fprintf(output, "      Efficiency | %16g +/- %-16g\n", W(efficiency));
  #undef W
}

static inline void
print_rcvXprim_side(FILE* output, const enum side side, struct rcvXprim* rXp)
{
  #define W(Name) rXp->Name[side].E, rXp->Name[side].SE
  fprintf(output, "                 |              [Incoming]                ");
  fprintf(output, "              [Absorbed]               \n");
  fprintf(output, "            Flux | %16g +/- %-16g | %16g +/- %-16g\n",
    W(in.flux), W(abs.flux));
  fprintf(output, "   Material loss | %16g +/- %-16g | %16g +/- %-16g\n",
    W(in.flux_mat_loss), W(abs.flux_mat_loss));
  fprintf(output, "Atmospheric loss | %16g +/- %-16g | %16g +/- %-16g\n",
    W(in.flux_atm_loss), W(abs.flux_atm_loss));
  fprintf(output, "No Material loss | %16g +/- %-16g | %16g +/- %-16g\n",
    W(in.flux_no_mat_loss), W(abs.flux_no_mat_loss));
  fprintf(output, "  No Atmos. loss | %16g +/- %-16g | %16g +/- %-16g\n",
    W(in.flux_no_atm_loss), W(abs.flux_no_atm_loss));
  #undef W
}

static inline void
print_simulation(FILE* output, struct simul* simul)
{
  size_t i;

  /* Global results */
  #define W(Name) simul->Name.E, simul->Name.SE
  fprintf(output, " Overall results (#Samples = %zu)\n", simul->nsamps);
  fprintf(output, "------------------------------------------------------------");
  fprintf(output, "------------------------------\n");
  fprintf(output, "  Potential flux | %16g +/- %-16g\n", W(potential_flux));
  fprintf(output, "   Absorbed flux | %16g +/- %-16g\n", W(absorbed_flux));
  fprintf(output, "   Cosine factor | %16g +/- %-16g\n", W(cos_factor));
  fprintf(output, "     Shadow loss | %16g +/- %-16g\n", W(shadow_loss));
  fprintf(output, "    Missing loss | %16g +/- %-16g\n", W(missing_loss));
  fprintf(output, "  Materials loss | %16g +/- %-16g\n", W(materials_loss));
  fprintf(output, "Atmospheric loss | %16g +/- %-16g\n", W(atmospheric_loss));
  fprintf(output, "\n");
  #undef W
  /* Per receivers results */
  FOR_EACH(i, 0, BUF_SZ(simul->rcvs)) {
    struct rcv* rcv = &BUF_AT(simul->rcvs, i);
    fprintf(output, " Receiver `%s' (Area = %g)\n", rcv->name, rcv->area);
    if(rcv->in.flux[FRONT].E >= 0) {
      fprintf(output, "-----------------------------------------------------[Front]");
      fprintf(output, "------------------------------\n");
      print_rcv_side(output, FRONT, rcv);
    }
    if(rcv->in.flux[BACK].E >= 0) {
      fprintf(output, "------------------------------------------------------[Back]");
      fprintf(output, "------------------------------\n");
      print_rcv_side(output, BACK, rcv);
    }
    fprintf(output, "\n");
  }
  /* Per primary results */
  FOR_EACH(i, 0, BUF_SZ(simul->prims)) {
    struct prim* prim = &BUF_AT(simul->prims, i);
    fprintf(output, " Primary `%s' (Area = %g; #Samples = %zu)\n",
      prim->name, prim->area, prim->nsamps);
    #define W(Name) prim->Name.E, prim->Name.SE
    fprintf(output, "------------------------------------------------------------");
    fprintf(output, "------------------------------\n");
    fprintf(output, "   Cosine factor | %16g +/- %-16g\n", W(cos_factor));
    fprintf(output, "     Shadow loss | %16g +/- %-16g\n", W(shadow_loss));
    fprintf(output, "\n");
    #undef W
  }
  /* Per receiverXprimary results */
  FOR_EACH(i, 0, BUF_SZ(simul->rcvXprims)) {
    struct rcvXprim* rXp = &BUF_AT(simul->rcvXprims, i);
    struct rcv* rcv = find_receiver_by_id(simul, rXp->rcv_id);
    struct prim* prim = find_primary_by_id(simul, rXp->prim_id);
    fprintf(output, " Receiver `%s' X Primary `%s'\n", rcv->name, prim->name);
    if(rXp->in.flux[FRONT].E >= 0) {
      fprintf(output, "-----------------------------------------------------[Front]");
      fprintf(output, "------------------------------\n");
      print_rcvXprim_side(output, FRONT, rXp);
    }
    if(rXp->in.flux[BACK].E >= 0) {
      fprintf(output, "------------------------------------------------------[Back]");
      fprintf(output, "------------------------------\n");
      print_rcvXprim_side(output, BACK, rXp);
    }
    fprintf(output, "\n");
  }
}

int
main(int argc, char** argv)
{
  char s[128];
  buf_char_T buf = BUF_NULL;
  FILE* input = stdin;
  char* line = NULL;

  if(argc > 1 && !(input = fopen(argv[1], "r"))) {
    fprintf(stderr, "Could not open the file `%s'.\n", argv[1]);
    return 1;
  }
  while((line = read_line(&buf, input))) {
    struct simul simul;
    FILE* output = NULL;

    if(strncmp(line, "#--- Sun direction:", 19)) continue;

    simul_init(&simul);

    CHK(!strncmp(line, "#--- Sun direction:", 19));
    CHK(sscanf(line+19, "%lf %lf (%*f %*f %*f)",
      &simul.azimuth, &simul.elevation)==2);
    CHK(snprintf(s, sizeof(s), "%g-%g-raw-results.txt",
      simul.azimuth, simul.elevation) < sizeof(s));
    read_simulation(&simul, input);

    printf("Writing `%s'\n", s);
    CHK(output = fopen(s, "w"));
    print_simulation(output, &simul);
    fclose(output);

    simul_release(&simul);
  }

  BUF_RELEASE(buf);
  if(input && input != stdin) fclose(input);
  return 0;
}

