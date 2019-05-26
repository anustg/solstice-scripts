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

struct mesh {
  BUF(double) coords;
  BUF(size_t) ids;
  BUF(size_t) entities;
  BUF(size_t) ncells;
  size_t voffset;
};
static const struct mesh MESH_NULL = {BUF_NULL, BUF_NULL, BUF_NULL, BUF_NULL, 0};

static inline void
mesh_release(struct mesh* msh)
{
  BUF_RELEASE(msh->coords);
  BUF_RELEASE(msh->ids);
  BUF_RELEASE(msh->entities);
  BUF_RELEASE(msh->ncells);
  memset(msh, 0, sizeof(struct mesh));
}

static inline void
mesh_write_vtk(FILE* output, struct mesh* msh)
{
  size_t i, n;

  fprintf(output, "# vtk DataFile Version 2.0\n");
  fprintf(output, "Test\n");
  fprintf(output, "ASCII\n");
  fprintf(output, "DATASET POLYDATA\n");

  n = BUF_SZ(msh->coords)/3;
  fprintf(output, "POINTS %zu double\n", n);
  FOR_EACH(i, 0, n) {
    fprintf(output, "%g %g %g\n",
      BUF_AT(msh->coords, i*3+0),
      BUF_AT(msh->coords, i*3+1),
      BUF_AT(msh->coords, i*3+2));
  }
  n = BUF_SZ(msh->ids)/3;
  fprintf(output, "POLYGONS %zu %zu\n", n, n*4);
  FOR_EACH(i, 0, n) {
    fprintf(output, "3 %zu %zu %zu\n",
      BUF_AT(msh->ids, i*3+0),
      BUF_AT(msh->ids, i*3+1),
      BUF_AT(msh->ids, i*3+2));
  }
}

static inline void
mesh_write_obj(FILE* output, struct mesh* msh)
{
  size_t i, n;

  n = BUF_SZ(msh->coords)/3;
  FOR_EACH(i, 0, n) {
    fprintf(output, "v %g %g %g\n",
      BUF_AT(msh->coords, i*3+0),
      BUF_AT(msh->coords, i*3+1),
      BUF_AT(msh->coords, i*3+2));
  }

  n = BUF_SZ(msh->ids)/3;
  FOR_EACH(i, 0, n) {
    fprintf(output, "f %zu %zu %zu\n",
      BUF_AT(msh->ids, i*3+0) + 1,
      BUF_AT(msh->ids, i*3+1) + 1,
      BUF_AT(msh->ids, i*3+2) + 1);
  }
}

static inline void
mesh_write_prim_data_vtk
  (FILE* output, const struct mesh* msh, struct simul* simul)
{
  struct prim* prim;
  struct rcvXprim* rcvXprim;
  size_t ircv;
  size_t iprim;
  size_t icell;
  size_t n;

  n = BUF_SZ(msh->ids) / 3;
  fprintf(output, "CELL_DATA %zu\n", n);

  fprintf(output, "FIELD PrimaryData %zu\n", 2 + BUF_SZ(simul->rcvs)*6);
  fprintf(output, "cos_factor 2 %zu double\n", n);
  FOR_EACH(iprim, 0, BUF_SZ(msh->entities)) {
    CHK(prim = find_primary_by_id(simul, BUF_AT(msh->entities, iprim)));
    FOR_EACH(icell, 0, BUF_AT(msh->ncells, iprim)) {
      fprintf(output, "%g %g\n", prim->cos_factor.E, prim->cos_factor.SE);
    }
  }

  fprintf(output, "shadow_loss 2 %zu double\n", n);
  FOR_EACH(iprim, 0, BUF_SZ(msh->entities)) {
    CHK(prim = find_primary_by_id(simul, BUF_AT(msh->entities, iprim)));
    FOR_EACH(icell, 0, BUF_AT(msh->ncells, iprim)) {
      fprintf(output, "%g %g\n", prim->shadow_loss.E, prim->shadow_loss.SE);
    }
  }

  #define WRITE(Side, Name) {                                                  \
    fprintf(output, "%s_"STR(Side)"_"STR(Name)" 2 %zu double\n", rcv->name, n);\
    FOR_EACH(iprim, 0, BUF_SZ(msh->entities)) {                                \
      CHK(rcvXprim = find_rcvXprim(simul,rcv->id,BUF_AT(msh->entities,iprim)));\
      FOR_EACH(icell, 0, BUF_AT(msh->ncells, iprim)) {                         \
        fprintf(output, "%g %g\n",                                             \
          rcvXprim->Name[Side].E,                                              \
          rcvXprim->Name[Side].SE);                                            \
      }                                                                        \
    }                                                                          \
  } (void)0
  FOR_EACH(ircv, 0, BUF_SZ(simul->rcvs)) {
    const struct rcv* rcv = &BUF_AT(simul->rcvs, ircv);
    WRITE(FRONT, in.flux);
    WRITE(FRONT, in.flux_mat_loss);
    WRITE(FRONT, in.flux_atm_loss);
    WRITE(BACK,  in.flux);
    WRITE(BACK,  in.flux_mat_loss);
    WRITE(BACK,  in.flux_atm_loss);
  }
  #undef WRITE
}

static inline void
mesh_write_rcv_data_vtk
  (FILE* output, const struct mesh* msh, struct simul* simul)
{
  struct rcv* rcv;
  size_t ircv;
  size_t icell;
  size_t n;

  n = BUF_SZ(msh->ids)/3;
  fprintf(output, "CELL_DATA %zu\n", n);
  fprintf(output, "FIELD PrimaryData 12\n");

  #define WRITE(Side, Name) {                                                  \
    fprintf(output, STR(Side)"_"STR(Name)" 2 %zu double\n", n);                \
    FOR_EACH(ircv, 0, BUF_SZ(msh->entities)) {                                 \
      CHK(rcv = find_receiver_by_id(simul, BUF_AT(msh->entities, ircv)));      \
      FOR_EACH(icell, 0, BUF_AT(msh->ncells, ircv)) {                          \
        fprintf(output, "%g %g\n", rcv->Name[Side].E, rcv->Name[Side].SE);     \
      }                                                                        \
    }                                                                          \
  } (void)0
  WRITE(FRONT, in.flux);
  WRITE(FRONT, in.flux_mat_loss);
  WRITE(FRONT, in.flux_atm_loss);
  WRITE(FRONT, efficiency);
  WRITE(BACK,  in.flux);
  WRITE(BACK,  in.flux_mat_loss);
  WRITE(BACK,  in.flux_atm_loss);
  WRITE(BACK,  efficiency);
  #undef WRITE

  #define WRITE_MAP(Flux, Side) {                                              \
    fprintf(output, STR(Side)"_"STR(Flux)"_map 2 %zu double\n", n);            \
    FOR_EACH(ircv, 0, BUF_SZ(msh->entities)) {                                 \
      CHK(rcv = find_receiver_by_id(simul, BUF_AT(msh->entities, ircv)));      \
      if(!BUF_SZ(rcv->map[Flux][Side])) {                                      \
        FOR_EACH(icell, 0, BUF_AT(msh->ncells,ircv)) fprintf(output,"-1 -1\n");\
      } else {                                                                 \
        FOR_EACH(icell, 0, BUF_AT(msh->ncells,ircv)) {                         \
          const struct mc* mc = &BUF_AT(rcv->map[Flux][Side], icell);          \
          fprintf(output, "%g %g\n", mc->E, mc->SE);                           \
        }                                                                      \
      }                                                                        \
    }                                                                          \
  } (void)0
  WRITE_MAP(INCOMING, FRONT);
  WRITE_MAP(INCOMING, BACK);
  WRITE_MAP(ABSORBED, FRONT);
  WRITE_MAP(ABSORBED, BACK);
  #undef WRITE_MAP
}

int
main(int argc, char** argv)
{
  buf_char_T buf = BUF_NULL;
  char* line = NULL;

  FILE* input;
  FILE* geom;
  FILE* output;

  if(argc < 3) {
    fprintf(stderr, "Usage: %s solstice-geometry solstice-simulation\n", argv[0]);
    return 1;
  }

  CHK(geom = fopen(argv[1], "r"));
  CHK(input = fopen(argv[2], "r"));

  CHK(read_line(&buf, geom)); /* Skip the sun direction of the geometry */
  while((line = read_line(&buf, input))) {
    struct simul simul;

    const struct rcv*  rcv  = NULL;
    const struct prim* prim = NULL;

    struct mesh msh_rcv  = MESH_NULL;
    struct mesh msh_prim = MESH_NULL;
    struct mesh msh_misc = MESH_NULL;

    char filename[128];
    char prefix[64];
    size_t ntris_grp = 0, nverts_grp = 0;
    size_t nverts_obj = 0, off_obj = 0, off_grp = 0;

    simul_init(&simul);

    CHK(!strncmp(line, "#--- Sun direction:", 19));
    CHK(sscanf(line+19, "%lf %lf (%*f %*f %*f)", &simul.azimuth, &simul.elevation)==2);
    CHK(snprintf(prefix, sizeof(prefix), "%g-%g-",
      simul.azimuth, simul.elevation) < sizeof(prefix));
    read_simulation(&simul, input);

    while((line = read_line(&buf, geom)) && strncmp(line, "#--- Sun", 8)) {
      if(!strncmp(line, "g ", 2)) {
        if(prim) {
          BUF_PUSH(msh_prim.ncells, ntris_grp);
          msh_prim.voffset += nverts_grp;
        }
        if(rcv) {
          BUF_PUSH(msh_rcv.ncells, ntris_grp);
          msh_rcv.voffset += nverts_grp;
        }
        if(!prim && !rcv) msh_misc.voffset += nverts_grp;
        prim = find_primary(&simul, line+2);
        rcv = find_receiver(&simul, line+2);
        if(prim) BUF_PUSH(msh_prim.entities, prim->id);
        if(rcv) BUF_PUSH(msh_rcv.entities, rcv->id);
        ntris_grp = 0;
        nverts_grp = 0;
        off_grp = 0;
        off_obj = nverts_obj;

      } else if(!strncmp(line, "v ", 2)) {
        double pos[3];
        CHK(sscanf(line+2, "%lf %lf %lf", pos+0, pos+1, pos+2) == 3);
        if(prim) {
          BUF_PUSH(msh_prim.coords, pos[0]);
          BUF_PUSH(msh_prim.coords, pos[1]);
          BUF_PUSH(msh_prim.coords, pos[2]);
        }
        if(rcv) {
          BUF_PUSH(msh_rcv.coords, pos[0]);
          BUF_PUSH(msh_rcv.coords, pos[1]);
          BUF_PUSH(msh_rcv.coords, pos[2]);
        }
        if(!prim && !rcv) {
          BUF_PUSH(msh_misc.coords, pos[0]);
          BUF_PUSH(msh_misc.coords, pos[1]);
          BUF_PUSH(msh_misc.coords, pos[2]);
        }
        ++nverts_grp;
        ++nverts_obj;
      } else if(!strncmp(line, "f ", 2)) {
        size_t tri[3];
        CHK(sscanf(line+2, "%zu %zu %zu", tri+0, tri+1, tri+2) == 3);
        if(prim) {
          BUF_PUSH(msh_prim.ids, tri[0]-1 + msh_prim.voffset + off_grp - off_obj);
          BUF_PUSH(msh_prim.ids, tri[1]-1 + msh_prim.voffset + off_grp - off_obj);
          BUF_PUSH(msh_prim.ids, tri[2]-1 + msh_prim.voffset + off_grp - off_obj);
        }
        if(rcv) {
          BUF_PUSH(msh_rcv.ids, tri[0]-1 + msh_rcv.voffset + off_grp - off_obj);
          BUF_PUSH(msh_rcv.ids, tri[1]-1 + msh_rcv.voffset + off_grp - off_obj);
          BUF_PUSH(msh_rcv.ids, tri[2]-1 + msh_rcv.voffset + off_grp - off_obj);
        }
        if(!prim && !rcv) {
          BUF_PUSH(msh_misc.ids, tri[0]-1 + msh_misc.voffset + off_grp - off_obj);
          BUF_PUSH(msh_misc.ids, tri[1]-1 + msh_misc.voffset + off_grp - off_obj);
          BUF_PUSH(msh_misc.ids, tri[2]-1 + msh_misc.voffset + off_grp - off_obj);
        }
        ++ntris_grp;
      } else if(!strcmp(line, "---")) {
        nverts_obj = 0;
        off_obj = 0;
        off_grp = nverts_grp;
      }
    }
    if(prim) BUF_PUSH(msh_prim.ncells, ntris_grp);
    if(rcv) BUF_PUSH(msh_rcv.ncells, ntris_grp);

    CHK(snprintf(filename, sizeof(filename),
      "%sprimaries.vtk", prefix) < sizeof(prefix));
    printf("Writing `%s'\n", filename);
    CHK(output = fopen(filename, "w"));
    mesh_write_vtk(output, &msh_prim);
    mesh_write_prim_data_vtk(output, &msh_prim, &simul);
    fclose(output);

    CHK(snprintf(filename, sizeof(filename),
      "%sreceivers.vtk", prefix) < sizeof(prefix));
    printf("Writing `%s'\n", filename);
    CHK(output = fopen(filename, "w"));
    mesh_write_vtk(output, &msh_rcv);
    mesh_write_rcv_data_vtk(output, &msh_rcv, &simul);
    fclose(output);

    CHK(snprintf(filename, sizeof(filename),
      "%smiscellaneous.obj", prefix) < sizeof(prefix));
    printf("Writing `%s'\n", filename);
    CHK(output = fopen(filename, "w"));
    mesh_write_obj(output, &msh_misc);
    fclose(output);

    mesh_release(&msh_prim);
    mesh_release(&msh_rcv);
    mesh_release(&msh_misc);
    simul_release(&simul);
  }

  fclose(geom);
  fclose(input);
  BUF_RELEASE(buf);
  return 0;
}

