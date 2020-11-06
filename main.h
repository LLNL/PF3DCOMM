
#define NBVMAX 6
#define NRVMAX 8
#define NPLANE 3

extern int  nodes, mp_p, mp_q, mp_r, mp_myp, mp_myq, mp_myr, mp_rank;
extern int  nxl, nyl, nzl;
extern real *arry, *brry, *zrry;
static real *inbuf= 0, *outbuf= 0, *z_inbuf= 0, *z_outbuf= 0;
static real *y_inbuf= 0, *y_outbuf= 0;
static real *x_inbuf= 0, *x_outbuf= 0;
static long bufsize=0;

extern MPI_Request  shiftRequest[6][2];

extern long bytes_moved;

extern void torows(real *checker, real *rows, int nx, int ny);
extern void fromrows(real *checker, real *rows, int nx, int ny);
extern void tocolumns(real *checker, real *cols, int nx, int ny);
extern void fromcolumns(real *checker, real *cols, int nx, int ny);
extern void torowsall(real *checker, real *rows, int nx, int ny);
extern void fromrowsall(real *rows, real *checker, int nx, int ny);
extern void tocolsall(real *checker, real *cols, int nx, int ny);
extern void fromcolsall(real *rows, real *checker, int nx, int ny);
extern void syncplane(real *checker, int nx, int ny, int nplane);
extern void syncplane_hydro(real *checker, int nx, int ny, int nz);
extern int storage_init(long numbuf);
