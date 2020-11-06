/*
 * g3transpose.c
 *
 * this file contains functions that transpose decomposed 2D arrays
 * from checkerboard to rows or columns and back using MPI.
 *
 */

#include "grid.h"

#include <stdio.h>
#include "main.h"

long bytes_moved= 0;


/*------------------------------------------------------------------------
 * build_grid
 *
 * define a 3-D Cartesian topology atop the hardware
 */

#define REORDER FALSE

/*
   These functions use four buffers. One has the original data and is used
   when the data is in checkerboard format and is called "checker".
   The second is used when the data is in either row or column format
   and is called "work". The third is called "outbuf" and holds outbound
   messages and the fourth is called "inbuf" and receives incoming messages.
   All 4 buffers are the same size.
  
   A transpose involves several shifts (mp_p-1 or mp_q-1, to be exact).
   Before the start of a transpose the data is copied from checker or work
   to outbuf. If converting from checker to work, bytes are copied from inbuf
   to work after each shift. If going in the opposite direction, bytes are copied
   from inbuf to checker after each shift.
   
   At the start of a transpose the bytes are copied from work or checker or
   work in such an order that a process copies bytes from the end of inbuf to
   their destination and then sends the initial part of inbuf out as the next
   message. When going from checker to work, rows or columns are rotated until the
   portions of the rows or columns needed by this process are at the end.
   The message is then passed to the process "before" this one, which then
   finds the bytes it needs at the end and so on. The pre-rotation means that
   the necessary bytes are always at the end of an incoming message and the
   bytes for the next (smaller) message are in contiguous memory.
   This is much simpler to work with than trying to grab bytes out of the middle
   of an incoming message and then performing a compaction before sending
   the next message out.

   The approach is somewhat different when going from work back to checker.
   When going from rows back to checker, the array is transposed from row index
   varying fastest to column index varying fastest. Each receiving process can then
   get the bytes it needs from the end of the message. The natural handling of
   this data results in checker having column index varying fastest.
   This is not the order expected by pf3d, but that doesn't matter - it is only
   in this order during an intermediate state of the 2D FFT. The final transpose
   back from columns back to checker involves first switching the row index varying
   fastest and then starting to send messages.

   IMPORTANT NOTE: This approach does not pass a minimal number of bytes - that
   would be done by sending the required bytes directly to every "recipient".
   It does minimize the number of bytes that must be sent in a "shift right"
   approach - at each stage only the bytes needed by "downwind" processes are
   sent out. This results in about half as many bytes being sent as in Bert Still's
   original "shift right" function where a full "square" from the checkerboard is
   passed at each step.
*/

/*------------------------------------------------------------------------
 * torows
 *
 * Transpose from a checkerboard decomposition to a decomposition by rows.
 * 
 */

void torows(real *checker, real *rows, int nx, int ny)
{
  MPI_Status is;
  long currlen, i, j, ip, curr_ip, rows_per_proc, size, chunk, base, ixlo;
  real *inptr, *outptr, *ptrtmp, *src, *dst;

  bytes_moved= 0;  /* count of bytes moved by this function */
  rows_per_proc= ny/g3.P;  /* number of rows per processor in the row decomp */
  /* Messages are sent to the process "before" this one in the
     x-decomposition because then the bytes that process needs will be at
     the end of the buffer it receives. */
  /* full size of a chunk owned by a process */
  size= currlen= nx*ny;
  /* number of elements in each message destined for this process */
  chunk= nx*rows_per_proc;
  /* position in the i-decomposition of the process from which data will next
     be received */
  curr_ip= g3.myP;
  /* Copy data from the input array to the "out" buffer. Rotate so that
     the rows destined for this process are last in the buffer.
     The modulus operation should be slow, so copy one row
     at a time.
  */
  base= ny-rows_per_proc;
  for(j= 0; j < ny; j++) {
    src= checker+j*nx;
    dst= outbuf+(j+base)%ny*nx;
    for(i= 0; i < nx; i++) {
      dst[i]= src[i];
    }
  }
  /* Copy the data already on this node into place */
  base= currlen-chunk;
  ixlo= curr_ip*nx;
  for(j= 0; j < rows_per_proc; j++) {
    for(i= 0; i < nx; i++) {
      rows[ixlo+i+j*nx]= outbuf[base+i+j*nx];
    }
  }
  inptr= inbuf;
  outptr= outbuf;
#if 0
  printf("On rank %d, will shift in X on %d processes\n", g3.me, g3.P);
#endif
  for(ip= 1; ip < g3.P; ip++) {
    currlen -= chunk;
    /* On each pass the message originates from the processor one further towards
       larger X than the last pass (modulo the number in the decomp). */
    curr_ip= (curr_ip+1)%g3.P;
    MPI_Sendrecv(inbuf, currlen, MPI_SINGLE, g3.xnbr[AFTER], SHIFT_TAG,
		 outbuf, currlen, MPI_SINGLE, g3.xnbr[BEFORE], SHIFT_TAG,
		 g3.xcom, &is);
    bytes_moved += currlen*sizeof(real);
    base= currlen-chunk;
    ixlo= curr_ip*nx;
    for(j= 0; j < rows_per_proc; j++) {
      for(i= 0; i < nx; i++) {
        rows[ixlo+i+j*nx]= inbuf[base+i+j*nx];
      }
    }
    /* swap buffers */
    ptrtmp= inptr;
    inptr= outptr;
    outptr= ptrtmp;
  }
}

/*------------------------------------------------------------------------
 * fromrows
 *
 * Transpose from a decomposition by rows to a checkerboard decomposition.
 * 
 */

void fromrows(real *checker, real *rows, int nx, int ny)
{
  MPI_Status is;
  long currlen, i, j, ip, curr_ip, rows_per_proc, size, chunk, base, ixlo, rowsize;
  long srcoff, dstoff, rowoff;
  real *inptr, *outptr, *ptrtmp, *src, *dst;

  bytes_moved= 0;  /* count of bytes moved by this function */
  rows_per_proc= ny/g3.P;  /* number of rows per processor in the row decomp */
  /* Messages are sent to the process "before" this one in the
     x-decomposition because then the bytes that process needs will be at
     the end of the buffer it receives. */
  /* full size of a chunk owned by a process */
  size= currlen= nx*ny;
  rowsize= nx*g3.P;
  /* number of elements in each message destined for this process */
  chunk= nx*rows_per_proc;
  /* position in the i-decomposition of the process from which data will next
     be received */
  curr_ip= g3.myP;
  /* Copy data from the rows array to the "out" buffer.
     Split rows into chunks of length nx and put the chunks destined
     for each processor together. "Rotate" the rows_per_proc chunks of
     length nx destined for this process to the end.
  */
  for(ip= 0; ip < g3.P; ip++) {
    srcoff= ip*nx;
    dstoff= (rowsize-nx+nx*(ip-g3.myP))%rowsize;
    /* gather all the data needed by process "ip" */
    for(j= 0; j < rows_per_proc; j++) {
      rowoff= j*rowsize;
      for(i= 0; i < nx; i++) {
        outbuf[i+dstoff+rowoff]= rows[i+srcoff+rowoff];
      }
    }
  }
  /* Copy the data already on this node into place */
  srcoff= currlen-chunk;
  dstoff= chunk*curr_ip;
  for(j= 0; j < rows_per_proc; j++) {
    rowoff= j*nx;
    for(i= 0; i < nx; i++) {
      checker[i+rowoff+dstoff]= outbuf[i+rowoff+srcoff];
    }
  }
  inptr= inbuf;
  outptr= outbuf;
#if 0
  if(g3.myP == 0) {
    printf("On rank 0, will shift in X on %d processes\n", g3.P);
  }
#endif
  for(ip= 1; ip < g3.P; ip++) {
    currlen -= chunk;
    /* On each pass the message originates from the processor one further towards
       larger X than the last pass (modulo the number in the decomp). */
    curr_ip= (curr_ip+1)%g3.P;
    MPI_Sendrecv(inbuf, currlen, MPI_SINGLE, g3.xnbr[AFTER], SHIFT_TAG,
		 outbuf, currlen, MPI_SINGLE, g3.xnbr[BEFORE], SHIFT_TAG,
		 g3.xcom, &is);
    bytes_moved += currlen*sizeof(real);
    srcoff= currlen-chunk;
    dstoff= chunk*curr_ip;
    for(j= 0; j < rows_per_proc; j++) {
      rowoff= j*nx;
      for(i= 0; i < nx; i++) {
        checker[i+rowoff+dstoff]= outbuf[i+rowoff+srcoff];
      }
    }
    /* swap buffers */
    ptrtmp= inptr;
    inptr= outptr;
    outptr= ptrtmp;
  }
}

/*------------------------------------------------------------------------
 * tocolumns
 *
 * Transpose from a checkerboard decomposition to a decomposition by columns.
 * 
 */
static int *requests, *statuses;

void tocolumns(real *checker, real *cols, int nx, int ny)
{
  MPI_Status is;
  long currlen, i, j, iq, curr_iq, cols_per_proc, size, chunk, base, iylo, in, jn;
  long srcoff, dstoff, coloff, colsize, ndx1, ndx2, ndx1p, ndx2p;
  real *inptr, *outptr, *ptrtmp, *src, *dst;

  bytes_moved= 0;  /* count of bytes moved by this function */
  cols_per_proc= nx/g3.Q;  /* number of cols per processor in the row decomp */
  /* Messages are sent to the process "before" this one in the
     y-decomposition because then the bytes that process needs will be at
     the end of the buffer it receives. */
  /* full size of a chunk owned by a process */
  size= currlen= nx*ny;
  colsize= ny*g3.Q;
  /* number of elements in each message destined for this process */
  chunk= ny*cols_per_proc;
  /* position in the j-decomposition of the process from which data will next
     be received */
  curr_iq= g3.myQ;
  /* Copy data from the input array to the "out" buffer. Transpose each tile
     and rotate so that the cols destined for this process are last in the buffer.
     The modulus operation should be slow, so copy one row at a time.
  */
  for(i= 0; i < nx; i++) {
    for(j= 0; j < ny; j++) {
      ndx1= j;
      ndx2= i;
      ndx1p= ndx1;
      ndx2p= (ndx2+cols_per_proc*(g3.Q-g3.myQ) )%ny;
      /* transpose and rotate at the same time */
      outbuf[ndx1p+ndx2p*ny]= checker[i+j*nx];
    }
  }
  /* Copy the data already on this node into place */
  srcoff= currlen-chunk;
  dstoff= chunk*curr_iq;
  for(i= 0; i < cols_per_proc; i++) {
    coloff= i*ny;
    for(j= 0; j < ny; j++) {
      cols[j+coloff+dstoff]= outbuf[j+coloff+srcoff];
    }
  }
  inptr= inbuf;
  outptr= outbuf;
#if 0
  if(g3.myQ == 0) {
    printf("On rank 0, will shift in X on %d processes\n", g3.Q);
  }
#endif
  for(iq= 1; iq < g3.Q; iq++) {
    currlen -= chunk;
    /* On each pass the message originates from the processor one further towards
       larger Y than the last pass (modulo the number in the decomp). */
    curr_iq= (curr_iq+1)%g3.Q;
    MPI_Sendrecv(inbuf, currlen, MPI_SINGLE, g3.ynbr[AFTER], SHIFT_TAG,
		 outbuf, currlen, MPI_SINGLE, g3.ynbr[BEFORE], SHIFT_TAG,
		 g3.ycom, &is);
    bytes_moved += currlen*sizeof(real);
    srcoff= currlen-chunk;
    dstoff= chunk*curr_iq;
    for(i= 0; i < cols_per_proc; i++) {
      coloff= i*ny;
      for(j= 0; j < ny; j++) {
        cols[j+coloff+dstoff]= outbuf[j+coloff+srcoff];
      }
    }
    /* swap buffers */
    ptrtmp= inptr;
    inptr= outptr;
    outptr= ptrtmp;
  }
}

/*------------------------------------------------------------------------
 * fromcolumns
 *
 * Transpose from a decomposition by columns to a checkerboard decomposition.
 * 
 */

void fromcolumns(real *checker, real *cols, int nx, int ny)
{
  MPI_Status is;
  long currlen, i, j, iq, curr_iq, cols_per_proc, size, chunk, iylo, in, jn;
  long srcoff, dstoff, coloff, colsize;
  real *inptr, *outptr, *ptrtmp, *src, *dst;

  bytes_moved= 0;  /* count of bytes moved by this function */
  /* Create message passing buffers if needed */
  cols_per_proc= nx/g3.Q;  /* number of cols per processor in the row decomp */
  /* Messages are sent to the process "before" this one in the
     y-decomposition because then the bytes that process needs will be at
     the end of the buffer it receives. */
  /* full size of a chunk owned by a process */
  size= currlen= nx*ny;
  colsize= ny*g3.Q;
  /* number of elements in each message destined for this process */
  chunk= ny*cols_per_proc;
  /* position in the j-decomposition of the process from which data will next
     be received */
  curr_iq= g3.myQ;
  /* Copy data from the columns array to the "out" buffer.
     Split columns into chunks of length nx and put the chunks destined
     for each processor together. "Rotate" the columns_per_proc chunks of
     length ny destined for this process to the end.
  */
  for(iq= 0; iq < g3.Q; iq++) {
    srcoff= iq*ny;
    dstoff= (colsize-ny+ny*(iq-g3.myQ))%colsize;
    /* gather all the data needed by process "iq" */
    for(i= 0; i < cols_per_proc; i++) {
      coloff= i*colsize;
      for(j= 0; j < nx; j++) {
        outbuf[j+dstoff+coloff]= cols[j+srcoff+coloff];
      }
    }
  }
  /* Copy the data already on this node into place */
  srcoff= currlen-chunk;
  dstoff= chunk*curr_iq;
  for(i= 0; i < cols_per_proc; i++) {
    coloff= i*ny;
    for(j= 0; j < nx; j++) {
      checker[j+coloff+dstoff]= outbuf[j+coloff+srcoff];
    }
  }
  inptr= inbuf;
  outptr= outbuf;
  for(iq= 1; iq < g3.Q; iq++) {
    currlen -= chunk;
    /* On each pass the message originates from the processor one further towards
       larger Y than the last pass (modulo the number in the decomp). */
    curr_iq= (curr_iq+1)%g3.Q;
    MPI_Sendrecv(inbuf, currlen, MPI_SINGLE, g3.ynbr[AFTER], SHIFT_TAG,
		 outbuf, currlen, MPI_SINGLE, g3.ynbr[BEFORE], SHIFT_TAG,
		 g3.ycom, &is);
    bytes_moved += currlen*sizeof(real);
#if 0
    if(g3.me == master) {
        printf("On message pass %d, g3.Q is %d and moved %ld bytes\n", iq, g3.Q, currlen*sizeof(real));
    }
#endif    
    srcoff= currlen-chunk;
    dstoff= chunk*curr_iq;
    for(i= 0; i < cols_per_proc; i++) {
      coloff= i*ny;
      for(j= 0; j < nx; j++) {
        checker[j+coloff+dstoff]= outbuf[j+coloff+srcoff];
      }
    }
    /* swap buffers */
    ptrtmp= inptr;
    inptr= outptr;
    outptr= ptrtmp;
  }
}

/*------------------------------------------------------------------------
 * torowsall
 *
 * Transpose from a checkerboard decomposition to a decomposition by rows.
 * Send messages directly to recipients instead of using a shift "right".
 * 
 */

void torowsall(real *checker, real *rows, int nx, int ny)
{
  MPI_Status is;
  long i, j, iz, ip, rows_per_proc, chunk, srcoff, dstoff, rowsize;
  long zones_xy= nxl*nyl, offset;

  /* Create message passing buffers if needed */
  /* number of rows per processor in the row decomp */
  rows_per_proc= ny/g3.P;
  /* number of elements to send to each other process */
  chunk= nx*rows_per_proc;
  rowsize= nx*g3.P;
  /*
    MPI_Alltoall is ideally suited for performing a transpose to rows.
    The input checkerboard array has the data in the proper order to be viewed as
    g3.P buffers of length nx*rows_per_proc.
    The alltoall command will copy the data already on this processor into
    the proper place as well as putting the messages from the other processors
    at their location in the buffer.
  */
  for(iz= 0; iz < nzl; iz++) {
    offset= zones_xy*iz;
    MPI_Alltoall(checker+offset, chunk, MPI_SINGLE, inbuf+offset, chunk, MPI_SINGLE, g3.xcom);
  }

  /*DOES THERE NEED TO BE A WAIT_ALL??? */

  /* count of bytes moved by this function */
  bytes_moved= nzl*chunk*g3.P*sizeof(real);
  /* all data is now here, but it needs to be put into row order */
  for(iz= 0; iz < nzl; iz++) {
     for(ip= 0; ip < g3.P; ip++) {
      for(j= 0; j < rows_per_proc; j++) {
        srcoff= j*nx+ip*chunk+zones_xy*iz;
        dstoff= j*rowsize+ip*nx+zones_xy*iz;
        for(i= 0; i < nx; i++) {
          rows[i+dstoff]= inbuf[i+srcoff];
        }
      }
    }
  }
}

/*------------------------------------------------------------------------
 * fromrowsall
 *
 * Transpose from a decomposition by rows to a checkerboard decomposition.
 * Use MPI_Alltoall.
 * 
 */

void fromrowsall(real *rows, real *checker, int nx, int ny)
{
  MPI_Status is;
  long i, j, iz, ip, rows_per_proc, chunk, srcoff, dstoff, rowsize;
  long zones_xy= nxl*nyl, offset;

  /* number of rows per processor in the row decomp */
  rows_per_proc= ny/g3.P;
  /* number of elements to send to each other process */
  chunk= nx*rows_per_proc;
  rowsize= nx*g3.P;
  /*
    MPI_Alltoall is well suited for performing a transpose from rows
    to checkerboard.
    Convert from full rows per processor to a linear array in which
    the data destined for each process is contiguous.
    At that point, it is in the order need for the Alltoall.
    The data in the receive buffer is in exactly checkerboard order
    so it is received directly into that array.
  */
  for(iz= 0; iz < nzl; iz++) {
    for(ip= 0; ip < g3.P; ip++) {
      for(j= 0; j < rows_per_proc; j++) {
        dstoff= j*nx+ip*chunk+zones_xy*iz;
        srcoff= j*rowsize+ip*nx+zones_xy*iz;
        for(i= 0; i < nx; i++) {
          inbuf[i+dstoff]= rows[i+srcoff];
        }
      }
    }
  }

  for(iz= 0; iz < nzl; iz++) {
    offset= zones_xy*iz;
    MPI_Alltoall(inbuf+offset, chunk, MPI_SINGLE, checker+offset, chunk, MPI_SINGLE, g3.xcom);
  }

  /* count of bytes moved by this function */
  bytes_moved= nzl*chunk*g3.P*sizeof(real);
}

/*------------------------------------------------------------------------
 * tocolsall
 *
 * Transpose from a checkerboard decomposition to a decomposition by columns.
 * Send messages directly to recipients instead of using a shift "up".
 * 
 */

void tocolsall(real *checker, real *cols, int nx, int ny)
{
  MPI_Status is;
  long i, j, iz, iq, cols_per_proc, chunk, srcoff, dstoff, colsize;
  long zones_xy= nxl*nyl, offset;

  /* Create message passing buffers if needed */
  /* number of cols per processor in the column decomp */
  cols_per_proc= nx/g3.Q;
  /* number of elements to send to each other process */
  chunk= ny*cols_per_proc;
  colsize= ny*g3.Q;
  for(iz= 0; iz < nzl; iz++) {
    for(j= 0; j < ny; j++) {
      for(i= 0; i < nx; i++) {
        outbuf[j+i*ny+zones_xy*iz]= checker[i+j*nx+zones_xy*iz];
      }
    }
  }
  /*
    MPI_Alltoall is ideally suited for performing a transpose to columns.
    Once the input checkerboard array has been transposed to column varying fastest
    order, the data is in the proper order to be viewed as
    g3.Q buffers of length ny*cols_per_proc.
    The alltoall command will copy the data already on this processor into
    the proper place as well as putting the messages from the other processors
    at their location in the buffer.
  */
  for(iz= 0; iz < nzl; iz++) {
    offset= zones_xy*iz;
    MPI_Alltoall(outbuf+offset, chunk, MPI_SINGLE, inbuf+offset, chunk, 
                 MPI_SINGLE, g3.ycom);
  }

  /* count of bytes moved by this function */
  bytes_moved= nzl*chunk*g3.Q*sizeof(real);
  /* all data is now here, but it needs to be put into full columns */
  for(iz= 0; iz < nzl; iz++) {
    for(iq= 0; iq < g3.Q; iq++) {
      for(j= 0; j < cols_per_proc; j++) {
        srcoff= j*ny+iq*chunk+zones_xy*iz;
        dstoff= j*colsize+iq*ny+zones_xy*iz;
        for(i= 0; i < ny; i++) {
          cols[i+dstoff]= inbuf [i+srcoff];
        }
      }
    }
  }
}



/*------------------------------------------------------------------------
 * fromcolsall
 *
 * Transpose froma decomposition by columns to a checkerboard decomposition.
 * Use MPI_Alltoall.
 * 
 */

void fromcolsall(real *cols, real *checker, int nx, int ny)
{
  MPI_Status is;
  long i, j, iz, iq, cols_per_proc, chunk, srcoff, dstoff, colsize;
  long zones_xy= nxl*nyl, offset;

  /* Create message passing buffers if needed */
  /* number of cols per processor in the column decomp */
  cols_per_proc= nx/g3.Q;
  /* number of elements to send to each other process */
  chunk= ny*cols_per_proc;
  colsize= ny*g3.Q;
  /* transpose the columns into a 1D array with the elements
     destined for each process being contiguous. */
  for(iz= 0; iz < nzl; iz++) {
    for(iq= 0; iq < g3.Q; iq++) {
      for(j= 0; j < cols_per_proc; j++) {
        dstoff= j*ny+iq*chunk+zones_xy*iz;
        srcoff= j*colsize+iq*ny+zones_xy*iz;
        for(i= 0; i < ny; i++) {
          inbuf[i+dstoff]= cols[i+srcoff];
        }
      }
    }
  }
  /*
    MPI_Alltoall is well suited for performing a transpose from columns
    to a checkerboard. After the above transpose, the data is in order
    for an Alltoall to deposit it in the destination array as a column
    major checkerboard array.
  */
  for(iz= 0; iz < nzl; iz++) {
    offset= zones_xy*iz;
    MPI_Alltoall(inbuf+offset, chunk, MPI_SINGLE, outbuf+offset, chunk, MPI_SINGLE, g3.ycom);
  }

  /* transpose from column major to row major. */
  for(j= 0; j < ny; j++) {
    for(i= 0; i < nx; i++) {
      checker[i+j*nx]= outbuf[j+i*ny];
    }
  }
  /* count of bytes moved by this function */
  bytes_moved= nzl*chunk*g3.Q*sizeof(real);
}



/*------------------------------------------------------------------------
 * syncplane
 *
 * send nplane xy boundary planes to the neighbor at greater
 * z and  nplane-1 planes to the neighbor at lower z.
 * Receive nplane planes from the neighbor at lower z
 * and receive nplane-1 planes from the neighbor at higher z.
 */
void syncplane(real *checker, int nx, int ny, int nplane)
{
  MPI_Request ireq;
  MPI_Status is;
  long i, j, k, ndx;
  int nzone;

  /* number of elements to send to the other process */
  nzone= nx*ny;
  /* In the 6th order advection package in pf3d, 
     three planes of float complex data are passed downstream
     and two planes are passed upstream.
  */
  ndx= 0;
  for(k= 0; k < 2*nplane; k++) {
    for(j= 0; j < ny; j++) {
      for(i= 0; i < nx; i++) {
        z_outbuf[ndx++]= checker[i+j*nx];
      }
    }
  }
  /* send data to downstream neighbor */
  if ( g3.myR < g3.R-1 ) {
    MPI_Isend(z_outbuf, 2*nplane*nzone, MPI_SINGLE, g3.znbr[AFTER], SBS_TAG, g3.zcom, &ireq);
  }
  /* receive data from upstream neighbor */
  if ( g3.myR > 0 ) {
    MPI_Recv(z_inbuf, 2*nplane*nzone, MPI_SINGLE, g3.znbr[BEFORE], SBS_TAG, g3.zcom, &is);
  }

  /* NOTE: z_outbuf can't safely be re-used until MPI_Isend completes. */
  if ( g3.myR < g3.R-1 ) {
    MPI_Wait(&ireq, &is);
  }
  ndx= 0;
  for(k= 0; k < 2*(nplane-1); k++) {
    for(j= 0; j < ny; j++) {
      for(i= 0; i < nx; i++) {
        z_outbuf[ndx++]= checker[i+j*nx];
      }
    }
  }

  /* send data to upstream neighbor */
  if ( g3.myR > 0 ) {
    MPI_Isend(z_outbuf, 2*(nplane-1)*nzone, MPI_SINGLE, g3.znbr[BEFORE], SBS_TAG, g3.zcom, &ireq);
  }
  /* receive data from downstream neighbor */
  if ( g3.myR < g3.R-1 ) {
    MPI_Recv(z_inbuf, 2*(nplane-1)*nzone, MPI_SINGLE, g3.znbr[AFTER], SBS_TAG, g3.zcom, &is);
  }

  if ( g3.myR > 0 ) {
    MPI_Wait(&ireq, &is);
  }

  /* count of bytes moved by this function */
  bytes_moved= (2*nplane+2*(nplane-1))*nzone*sizeof(real);
}


/*------------------------------------------------------------------------
 * syncplane_hydro
 *
 * pass boundary planes to all 6 neighbors as is done in the hydro package
 */
void syncplane_hydro(real *checker, int nx, int ny, int nz)
{
  MPI_Status is;
  long i, j, k, ndx;
  int nzone;

  bytes_moved= 0;


  /* XY Planes */
  /* number of elements to send to the other process */
  nzone= nx*ny;
  /* Pass one xy plane of data to the "after" neighbor in z
     and receive from the "before" neighbor in z.
  */
  ndx= 0;
  for(j= 0; j < ny; j++) {
      for(i= 0; i < nx; i++) {
        z_outbuf[ndx++]= checker[i+j*nx];
      }
  }
  /* send data to downstream neighbor */
  if ( g3.myR < g3.R-1 ) {
    MPI_Send(z_outbuf, nzone, MPI_SINGLE, g3.znbr[AFTER], SBS_TAG, g3.zcom);
  }
  /* receive data from upstream neighbor */
  if ( g3.myR > 0 ) {
    MPI_Recv(z_inbuf, nzone, MPI_SINGLE, g3.znbr[BEFORE], SBS_TAG, g3.zcom, &is);
  }
  /* count of bytes moved by this function */
  bytes_moved += nzone*sizeof(real);

  /* Pass the same xy plane  to the "before" neighbor in z and
     receive from the "after" neighbor in z.
  */
  /* send data to upstream neighbor */
  if ( g3.myR > 0 ) {
    MPI_Send(z_outbuf, nzone, MPI_SINGLE, g3.znbr[BEFORE], SBS_TAG, g3.zcom);
  }
  /* receive data from downstream neighbor */
  if ( g3.myR < g3.R-1 ) {
    MPI_Recv(z_inbuf, nzone, MPI_SINGLE, g3.znbr[AFTER], SBS_TAG, g3.zcom, &is);
  }
  /* count of bytes moved by this function */
  bytes_moved += nzone*sizeof(real);


  /* XZ Planes */
  /* number of elements to send to the other process */
  nzone= nx*nz;
  /* Pass one xz plane of data to the "after" neighbor in y
     and receive from the "before" neighbor in y.
  */
  ndx= 0;
  for(j= 0; j < nz; j++) {
      for(i= 0; i < nx; i++) {
        y_outbuf[ndx++]= checker[i+j*nx];
      }
  }
  /* send data to downstream neighbor */
  if ( g3.myQ < g3.Q-1 ) {
    MPI_Send(y_outbuf, nzone, MPI_SINGLE, g3.ynbr[AFTER], SBS_TAG, g3.ycom);
  }
  /* receive data from upstream neighbor */
  if ( g3.myQ > 0 ) {
    MPI_Recv(y_inbuf, nzone, MPI_SINGLE, g3.ynbr[BEFORE], SBS_TAG, g3.ycom, &is);
  }
  /* count of bytes moved by this function */
  bytes_moved += nzone*sizeof(real);

  /* Pass the same xz plane  to the "before" neighbor in y and
     receive from the "after" neighbor in y.
  */
  /* send data to upstream neighbor */
  if ( g3.myQ > 0 ) {
    MPI_Send(y_outbuf, nzone, MPI_SINGLE, g3.ynbr[BEFORE], SBS_TAG, g3.ycom);
  }
  /* receive data from downstream neighbor */
  if ( g3.myQ < g3.Q-1 ) {
    MPI_Recv(y_inbuf, nzone, MPI_SINGLE, g3.ynbr[AFTER], SBS_TAG, g3.ycom, &is);
  }
  /* count of bytes moved by this function */
  bytes_moved += nzone*sizeof(real);


  /* YZ Planes */
  /* number of elements to send to the other process */
  nzone= ny*nz;
  /* Pass one yz plane of data to the "after" neighbor in x
     and receive from the "before" neighbor in x.
  */
  ndx= 0;
  for(j= 0; j < nz; j++) {
      for(i= 0; i < ny; i++) {
        x_outbuf[ndx++]= checker[i+j*ny];
      }
  }
  /* send data to downstream neighbor */
  if ( g3.myP < g3.P-1 ) {
    MPI_Send(x_outbuf, nzone, MPI_SINGLE, g3.xnbr[AFTER], SBS_TAG, g3.xcom);
  }
  /* receive data from upstream neighbor */
  if ( g3.myP > 0 ) {
    MPI_Recv(x_inbuf, nzone, MPI_SINGLE, g3.xnbr[BEFORE], SBS_TAG, g3.xcom, &is);
  }
  /* count of bytes moved by this function */
  bytes_moved += nzone*sizeof(real);

  /* Pass the same yz plane  to the "before" neighbor in x and
     receive from the "after" neighbor in x.
  */
  /* send data to upstream neighbor */
  if ( g3.myP > 0 ) {
    MPI_Send(x_outbuf, nzone, MPI_SINGLE, g3.xnbr[BEFORE], SBS_TAG, g3.xcom);
  }
  /* receive data from downstream neighbor */
  if ( g3.myP < g3.P-1 ) {
    MPI_Recv(x_inbuf, nzone, MPI_SINGLE, g3.xnbr[AFTER], SBS_TAG, g3.xcom, &is);
  }
  /* count of bytes moved by this function */
  bytes_moved += nzone*sizeof(real);
}



int storage_init(long numbuf)
{
  extern int  nxl, nyl, nzl;
  long bytalloc= 0, num_xyz, num_xy, num_xz, num_yz;
  long len_xyz, len_xy, len_xz, len_yz, byt2d, byt3d, byttot;
  
  /* Allocate all scratch space needed for this run.
     This function can be called to re-allocate storage.
  */
  num_xyz= 4;
  num_xy= 3;
  num_xz= 2;
  num_yz= 2;
  len_xyz= 2*nxl*nyl*nzl;
  len_xy= nxl*nyl;
  len_xz= nxl*nzl;
  len_yz= nyl*nzl;
  byt2d= (num_xy*len_xy+num_xz*len_xz +num_yz*len_yz)*NPLANE*sizeof(real)*2;
  byt3d= num_xyz*numbuf*sizeof(real);
  byttot= byt2d+byt3d;
  if(g3.me == master) printf("Expected number of bytes allocated is %ld\n", 
			byttot);
  if(bufsize < numbuf) {
    if(arry) free(arry);
    arry = (real *) calloc( numbuf, sizeof(real) );
    if(!arry) {
      bufsize= 0;
      return -1;
    }
    bytalloc += numbuf*sizeof(real);
    if(brry) free(brry);
    brry = (real *) calloc( numbuf, sizeof(real) );
    if(!brry) {
      bufsize= 0;
      return -1;
    }
    bytalloc += numbuf*sizeof(real);

    if(zrry) free(zrry);
    zrry = (real *) calloc( 2*nxl*nyl*NPLANE, sizeof(real) );
    if(!zrry) {
      bufsize= 0;
      return -1;
    }
    bytalloc += 2*nxl*nyl*NPLANE*sizeof(real);

    if(inbuf) free(inbuf);
    inbuf= (real *) malloc(numbuf*sizeof(real));
    if(!inbuf) {
      bufsize= 0;
      return -2;
    }
    bytalloc += numbuf*sizeof(real);
    if(outbuf) free(outbuf);
    outbuf= (real *) malloc(numbuf*sizeof(real));
    if(!outbuf) {
      bufsize= 0;
      return -2;
    }
    bytalloc += numbuf*sizeof(real);

    if(z_inbuf) free(z_inbuf);
    z_inbuf= (real *) malloc(2*nxl*nyl*NPLANE*sizeof(real));
    if(!z_inbuf) {
      bufsize= 0;
      return -2;
    }
    bytalloc += 2*nxl*nyl*NPLANE*sizeof(real);
    if(z_outbuf) free(z_outbuf);
    z_outbuf= (real *) malloc(2*nxl*nyl*NPLANE*sizeof(real));
    if(!z_outbuf) {
      bufsize= 0;
      return -2;
    }
    bytalloc += 2*nxl*nyl*NPLANE*sizeof(real);
    if(y_inbuf) free(y_inbuf);
    y_inbuf= (real *) malloc(2*nxl*nzl*sizeof(real));
    if(!y_inbuf) {
      bufsize= 0;
      return -2;
    }
    bytalloc += 2*nxl*nzl*NPLANE*sizeof(real);
    if(y_outbuf) free(y_outbuf);
    y_outbuf= (real *) malloc(2*nxl*nzl*sizeof(real));
    if(!y_outbuf) {
      bufsize= 0;
      return -2;
    }
    bytalloc += 2*nxl*nzl*NPLANE*sizeof(real);
    if(x_inbuf) free(x_inbuf);
    x_inbuf= (real *) malloc(2*nyl*nzl*sizeof(real));
    if(!x_inbuf) {
      bufsize= 0;
      return -2;
    }
    bytalloc += 2*nyl*nzl*NPLANE*sizeof(real);
    if(x_outbuf) free(x_outbuf);
    x_outbuf= (real *) malloc(2*nyl*nzl*sizeof(real));
    if(!x_outbuf) {
      bufsize= 0;
      return -2;
    }
    bytalloc += 2*nyl*nzl*NPLANE*sizeof(real);

    bufsize= numbuf;
    if(g3.me == master) printf("*** Number of bytes allocated is %ld\n", bytalloc);
  }
  return 0;
}
