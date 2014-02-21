#define INDEX2(xsize,x,y) x + (y)*xsize
#define INDEX3(xsize,ysize,x,y,z) x + xsize*(y + ysize*z)
#define INDEX4(xsize,ysize,zsize,x,y,z,i) x + xsize*(y + ysize*(z + zsize*i))
#define INDEX5(xsize,ysize,zsize,isize,x,y,z,i,j) x + xsize*(y + ysize*(z + zsize*(i + isize*(j))))
#define NDIM 3
#define NGLLX 5
#define NGLL2 25
#define NGLL3 125
#define NGLL3_PADDED 128
#define N_SLS 3
#define IREGION_CRUST_MANTLE 1
#define IREGION_INNER_CORE 3
#define IFLAG_IN_FICTITIOUS_CUBE 11
#define R_EARTH_KM 6371.0f
#define COLORING_MIN_NSPEC_INNER_CORE 1000
#define COLORING_MIN_NSPEC_OUTER_CORE 1000
#define BLOCKSIZE_TRANSFER 256
__global__ void write_seismograms_transfer_from_device_kernel(const int * number_receiver_global, const int * ispec_selected_rec, const int * ibool, float * station_seismo_field, const float * desired_field, const int nrec_local){
  int blockID;
  blockID = blockIdx.x + (blockIdx.y) * (gridDim.x);
  if(blockID < nrec_local){
    int irec;
    int ispec;
    int iglob;
    irec = number_receiver_global[blockID - 0] - (1);
    ispec = ispec_selected_rec[irec - 0] - (1);
    iglob = ibool[threadIdx.x + (NGLL3) * (ispec) - 0] - (1);
    station_seismo_field[((NGLL3) * (3)) * (blockID) + (threadIdx.x) * (3) + 0 - 0] = desired_field[(iglob) * (3) + 0 - 0];
    station_seismo_field[((NGLL3) * (3)) * (blockID) + (threadIdx.x) * (3) + 1 - 0] = desired_field[(iglob) * (3) + 1 - 0];
    station_seismo_field[((NGLL3) * (3)) * (blockID) + (threadIdx.x) * (3) + 2 - 0] = desired_field[(iglob) * (3) + 2 - 0];
  }
}
