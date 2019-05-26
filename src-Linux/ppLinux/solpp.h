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

#define _POSIX_C_SOURCE 200809L /* Support of strdup */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef SOLPP_H
#define SOLPP_H

enum side { FRONT, BACK };
enum flux_density { ABSORBED, INCOMING };
struct mc { double E/* Expected value */, SE/* Standard Error */; };

/*******************************************************************************
 * Helper macros
 ******************************************************************************/
#define STR__(X) #X
#define STR(X) STR__(X)

#define FOR_EACH(Id, Start, End) for((Id)=(Start); (Id)<(End); ++(Id))

#define CHK(Cond) if(!(Cond)) \
  {fprintf(stderr, "error:%s:%d\n", __FILE__, __LINE__); abort();} (void)0

/*******************************************************************************
 * Dynamic buffer
 ******************************************************************************/
#define BUF(Type) struct { Type* mem; size_t ca; size_t sz; }
#define BUF_NULL {NULL, 0, 0}
#define BUF_RELEASE(B) free((B).mem)
#define BUF_RESERVE(B, Sz)                                                     \
  if((Sz)>(B).ca) {                                                            \
    CHK((B).mem = realloc((B).mem, sizeof(*(B).mem)*((B).ca=(Sz))));           \
  } (void)0
#define BUF_RESIZE(B, Sz) {BUF_RESERVE((B), Sz); (B).sz = Sz;} (void)0
#define BUF_PUSH(B, E) {                                                       \
  if((B).sz >= (B).ca) BUF_RESERVE((B), ((B).ca?(B).ca*2:32));                 \
  (B).mem[(B).sz++]=(E);                                                       \
} (void)0
#define BUF_SZ(B) (B).sz
#define BUF_MEM(B) (B).mem
#define BUF_AT(B, I) (B).mem[I]

/*******************************************************************************
 * Helper function
 ******************************************************************************/
typedef BUF(char) buf_char_T;

static inline char*
read_line(buf_char_T* b, FILE* stream)
{
  char* c;

  if(!BUF_SZ(*b)) BUF_RESIZE(*b, 32);
  if(!fgets(BUF_MEM(*b), (int)BUF_SZ(*b), stream)) return NULL;

  /* Ensure that the whole line is read */
  while(!strrchr(BUF_MEM(*b), '\n') && !feof(stream)) {
    BUF_RESIZE(*b, BUF_SZ(*b)+32);
    CHK(fgets
      (BUF_MEM(*b) + strlen(BUF_MEM(*b)),
       (int)(BUF_SZ(*b) - strlen(BUF_MEM(*b))), stream));
  }

  /* Remove the carriage return */
  if((c = strrchr(BUF_MEM(*b), '\n'))) *c = '\0';

  return BUF_MEM(*b);
}

/*******************************************************************************
 * Per receiver double sided estimations
 ******************************************************************************/
struct flux {
  struct mc flux[2];
  struct mc flux_no_mat_loss[2];
  struct mc flux_no_atm_loss[2];
  struct mc flux_mat_loss[2];
  struct mc flux_atm_loss[2];
};

struct rcv {
  struct flux in;
  struct flux abs;
  struct mc efficiency[2];
  BUF(struct mc) map[2][2]; /* Absorbed/Incoming, Front/Back */
  char* name;
  size_t id;
  double area;
};

static inline void rcv_init(struct rcv* r) {memset(r, 0, sizeof(*r));}
static inline void rcv_release(struct rcv* r)
{
  free(r->name);
  BUF_RELEASE(r->map[ABSORBED][FRONT]);
  BUF_RELEASE(r->map[ABSORBED][BACK]);
  BUF_RELEASE(r->map[INCOMING][FRONT]);
  BUF_RELEASE(r->map[INCOMING][BACK]);
}

/*******************************************************************************
 * Per primary estimations
 ******************************************************************************/
struct prim {
  struct mc cos_factor;
  struct mc shadow_loss;
  char* name;
  size_t id;
  size_t nsamps;
  double area;
};
static inline void prim_init(struct prim* p) {memset(p, 0, sizeof(*p));}
static inline void prim_release(struct prim* p) {free(p->name);}

/*******************************************************************************
 * Per receiver X primary estimations
 ******************************************************************************/
struct rcvXprim {
  struct flux in;
  struct flux abs;
  size_t rcv_id;
  size_t prim_id;
};
static inline void rcvXprim_init(struct rcvXprim* r) {memset(r, 0, sizeof(*r));}
static inline void rcvXprim_release(struct rcvXprim* r) { /*Do nothing*/ }

/*******************************************************************************
 * Overall estimations of a simulation
 ******************************************************************************/
struct simul {
  struct mc potential_flux;
  struct mc absorbed_flux;
  struct mc cos_factor;
  struct mc shadow_loss;
  struct mc missing_loss;
  struct mc materials_loss;
  struct mc atmospheric_loss;

  BUF(struct rcv) rcvs;
  BUF(struct prim) prims;
  BUF(struct rcvXprim) rcvXprims;

  double azimuth;
  double elevation;
  size_t nsamps;
};

static inline void simul_init(struct simul* s) {memset(s, 0, sizeof(*s));}

static inline void
simul_release(struct simul* s)
{
  size_t i;
  FOR_EACH(i,0,BUF_SZ(s->rcvs)) rcv_release(&BUF_AT(s->rcvs, i));
  FOR_EACH(i,0,BUF_SZ(s->prims)) prim_release(&BUF_AT(s->prims, i));
  FOR_EACH(i,0,BUF_SZ(s->rcvXprims)) rcvXprim_release(&BUF_AT(s->rcvXprims, i));
  BUF_RELEASE(s->rcvs);
  BUF_RELEASE(s->prims);
  BUF_RELEASE(s->rcvXprims);
}

/*******************************************************************************
 * Look for entities
 ******************************************************************************/
static inline struct prim*
find_primary(struct simul* s, const char* name)
{
  size_t i;
  FOR_EACH(i, 0, BUF_SZ(s->prims))
    if(!strcmp(BUF_AT(s->prims, i).name, name)) return &BUF_AT(s->prims, i);
  return NULL;
}

static inline struct prim*
find_primary_by_id(struct simul* s, const size_t id)
{
  size_t i;
  FOR_EACH(i, 0, BUF_SZ(s->prims))
    if(BUF_AT(s->prims, i).id == id) return &BUF_AT(s->prims, i);
  return NULL;
}

static inline struct rcv*
find_receiver(struct simul* s, const char* name)
{
  size_t i;
  FOR_EACH(i, 0, BUF_SZ(s->rcvs))
    if(!strcmp(BUF_AT(s->rcvs, i).name, name)) return &BUF_AT(s->rcvs, i);
  return NULL;
}

static inline struct rcv*
find_receiver_by_id(struct simul* s, const size_t id)
{
  size_t i;
  FOR_EACH(i, 0, BUF_SZ(s->rcvs))
    if(BUF_AT(s->rcvs, i).id == id) return &BUF_AT(s->rcvs, i);
  return NULL;
}

static inline struct rcvXprim*
find_rcvXprim(struct simul* s, const size_t rcv_id, const size_t prim_id)
{
  size_t i;
  FOR_EACH(i, 0, BUF_SZ(s->rcvXprims)) {
    struct rcvXprim* rXp = &BUF_AT(s->rcvXprims, i);
    if(rXp->rcv_id == rcv_id && rXp->prim_id == prim_id) return rXp;
  }
  return NULL;
}

/*******************************************************************************
 * Read simulation data from a solstice-output
 ******************************************************************************/
static inline void
read_receiver_map_side_data(struct rcv* rcv, const size_t n, FILE* input)
{
  buf_char_T buf = BUF_NULL;
  char* line = NULL;
  size_t i;
  enum side side;
  enum flux_density flux;
  const char* str;

  CHK(line = read_line(&buf, input));
  str = line + 8;
       if(!strncmp(str, "Front_faces", 11)) { side = FRONT; str += 12; }
  else if(!strncmp(str, "Back_faces",  10)) { side = BACK;  str += 11; }
  else { fprintf(stderr, "Unexpected side name\n"); abort(); }

       if(!strncmp(str, "Incoming_flux", 13)) { flux = INCOMING; }
  else if(!strncmp(str, "Absorbed_flux", 13)) { flux = ABSORBED; }
  else { fprintf(stderr, "Unexpected flux name\n"); abort(); }

  CHK(read_line(&buf, input)); /* Discard "LOOKUP_TABLE default" line */

  BUF_RESIZE(rcv->map[flux][side], n);
  FOR_EACH(i, 0, n) {
    struct mc* mc = &BUF_AT(rcv->map[flux][side], i);
    CHK(line = read_line(&buf, input));
    CHK(sscanf(line, "%lf %lf", &mc->E, &mc->SE) == 2);
  }

  BUF_RELEASE(buf);
}

static inline void
read_receiver_map(struct simul* simul, FILE* input)
{
  struct rcv* rcv;
  buf_char_T buf = BUF_NULL;
  char* line = NULL;
  size_t i, n;
  long fp;

  CHK(line = read_line(&buf, input));
  CHK(rcv = find_receiver(simul, line));

  /* Skip header */
  CHK(read_line(&buf, input));
  CHK(read_line(&buf, input));
  /* Skip vertices */
  CHK(line = read_line(&buf, input));
  CHK(sscanf(line, "POINTS  %zu float", &n) == 1);
  FOR_EACH(i, 0, n) CHK(read_line(&buf, input));
  /* Skip polygons */
  CHK(line = read_line(&buf, input));
  CHK(sscanf(line, "POLYGONS %zu %*u", &n) == 1);
  FOR_EACH(i, 0, n) CHK(read_line(&buf, input));
  /* Read the map data of one side */
  CHK(line = read_line(&buf, input));
  CHK(sscanf(line, "CELL_DATA %zu", &n) == 1);
  /* Read map data */
  do {
    read_receiver_map_side_data(rcv, n, input);
    fp = ftell(input);
    line = read_line(&buf, input);
    fseek(input, fp, SEEK_SET);
  } while(line && !strncmp(line, "SCALARS", 7));

  BUF_RELEASE(buf);
}

static inline void
read_simulation(struct simul* simul, FILE* input)
{
  buf_char_T buf = BUF_NULL;
  char* line = NULL;
  char* tk = NULL;
  size_t nrcvs, nprims;
  size_t i;

  /* Counters */
  CHK(line = read_line(&buf, input));
  CHK(sscanf(line, "%*u %zu %zu %zu %*u", &nrcvs, &nprims, &simul->nsamps)==3);

  /* Global results */
  #define READ(Name) {                                                         \
    CHK(line = read_line(&buf, input));                                        \
    CHK(sscanf(line, "%lf %lf", &simul->Name.E, &simul->Name.SE) == 2);        \
  } (void)0
  READ(potential_flux);
  READ(absorbed_flux);
  READ(cos_factor);
  READ(shadow_loss);
  READ(missing_loss);
  READ(materials_loss);
  READ(atmospheric_loss);
  #undef READ

  /* Read per receiver results */
  BUF_RESIZE(simul->rcvs, nrcvs);
  FOR_EACH(i, 0, nrcvs) {
    struct rcv* rcv = &BUF_AT(simul->rcvs, i);
    rcv_init(rcv);

    CHK(line = read_line(&buf, input));
    CHK(tk = strtok(line, " \t"));
    CHK(rcv->name = strdup(tk));

    CHK(tk = strtok(NULL, ""));
    #define GET(Side, Name) &rcv->Name[Side].E, &rcv->Name[Side].SE
    CHK(sscanf
      (tk,
       "%zu %lf "
       "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
       "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
       "%lf %lf "
       "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
       "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
       "%lf %lf",
       &rcv->id, &rcv->area,
       GET(FRONT, in.flux),
       GET(FRONT, in.flux_no_mat_loss),
       GET(FRONT, in.flux_no_atm_loss),
       GET(FRONT, in.flux_mat_loss),
       GET(FRONT, in.flux_atm_loss),
       GET(FRONT, abs.flux),
       GET(FRONT, abs.flux_no_mat_loss),
       GET(FRONT, abs.flux_no_atm_loss),
       GET(FRONT, abs.flux_mat_loss),
       GET(FRONT, abs.flux_atm_loss),
       GET(FRONT, efficiency),
       GET(BACK, in.flux),
       GET(BACK, in.flux_no_mat_loss),
       GET(BACK, in.flux_no_atm_loss),
       GET(BACK, in.flux_mat_loss),
       GET(BACK, in.flux_atm_loss),
       GET(BACK, abs.flux),
       GET(BACK, abs.flux_no_mat_loss),
       GET(BACK, abs.flux_no_atm_loss),
       GET(BACK, abs.flux_mat_loss),
       GET(BACK, abs.flux_atm_loss),
       GET(BACK, efficiency)) == 46);
    #undef GET
  }

  /* Read per primary results */
  BUF_RESIZE(simul->prims, nprims);
  FOR_EACH(i, 0, nprims) {
    struct prim* prim = &BUF_AT(simul->prims, i);
    prim_init(prim);

    CHK(line = read_line(&buf, input));
    CHK(tk = strtok(line, " \t"));
    CHK(prim->name = strdup(tk));

    CHK(tk = strtok(NULL, ""));
    CHK(sscanf(tk, "%zu %lf %zu %lf %lf %lf %lf",
      &prim->id, &prim->area, &prim->nsamps,
      &prim->cos_factor.E, &prim->cos_factor.SE,
      &prim->shadow_loss.E, &prim->shadow_loss.SE) == 7);
  }

  /* Per receiverXprimary results */
  BUF_RESIZE(simul->rcvXprims, nprims*nrcvs);
  FOR_EACH(i, 0, nprims*nrcvs) {
    struct rcvXprim* rcvXprim = &BUF_AT(simul->rcvXprims, i);
    rcvXprim_init(rcvXprim);

    CHK(line = read_line(&buf, input));
    #define GET(Side, Name) &rcvXprim->Name[Side].E, &rcvXprim->Name[Side].SE
    CHK(sscanf
      (line,
       "%zu %zu "
       "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
       "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
       "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
       "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
       &rcvXprim->rcv_id, &rcvXprim->prim_id,
       GET(FRONT, in.flux),
       GET(FRONT, in.flux_no_mat_loss),
       GET(FRONT, in.flux_no_atm_loss),
       GET(FRONT, in.flux_mat_loss),
       GET(FRONT, in.flux_atm_loss),
       GET(FRONT, abs.flux),
       GET(FRONT, abs.flux_no_mat_loss),
       GET(FRONT, abs.flux_no_atm_loss),
       GET(FRONT, abs.flux_mat_loss),
       GET(FRONT, abs.flux_atm_loss),
       GET(BACK, in.flux),
       GET(BACK, in.flux_no_mat_loss),
       GET(BACK, in.flux_no_atm_loss),
       GET(BACK, in.flux_mat_loss),
       GET(BACK, in.flux_atm_loss),
       GET(BACK, abs.flux),
       GET(BACK, abs.flux_no_mat_loss),
       GET(BACK, abs.flux_no_atm_loss),
       GET(BACK, abs.flux_mat_loss),
       GET(BACK, abs.flux_atm_loss)) == 42);
    #undef GET
  }

  /* Read receiver maps */
  for(;;) {
    const long fp = ftell(input);
    line = read_line(&buf, input);

    if(line && !strncmp(line, "# vtk", 5)) {
      read_receiver_map(simul, input);
    } else {
      fseek(input, fp, SEEK_SET);
      break;
    }
  }
  BUF_RELEASE(buf);
}

#endif /* SOLPP_H */

