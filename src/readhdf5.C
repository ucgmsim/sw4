//  SW4 LICENSE
// # ----------------------------------------------------------------------
// # SW4 - Seismic Waves, 4th order
// # ----------------------------------------------------------------------
// # Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// # Produced at the Lawrence Livermore National Laboratory. 
// # 
// # Written by:
// # N. Anders Petersson (petersson1@llnl.gov)
// # Bjorn Sjogreen      (sjogreen2@llnl.gov)
// # 
// # LLNL-CODE-643337 
// # 
// # All rights reserved. 
// # 
// # This file is part of SW4, Version: 1.0
// # 
// # Please also read LICENCE.txt, which contains "Our Notice and GNU General Public License"
// # 
// # This program is free software; you can redistribute it and/or modify
// # it under the terms of the GNU General Public License (as published by
// # the Free Software Foundation) version 2, dated June 1991. 
// # 
// # This program is distributed in the hope that it will be useful, but
// # WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
// # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
// # conditions of the GNU General Public License for more details. 
// # 
// # You should have received a copy of the GNU General Public License
// # along with this program; if not, write to the Free Software
// # Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA 

#ifndef READHDF5_C
#define READHDF5_C

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <unistd.h>
#include <ctime>

#include "Require.h"
#include "EW.h"
#include "GridGenerator.h"
#include "TimeSeries.h"
#include "readhdf5.h"

#ifdef USE_HDF5
#include "hdf5.h"

struct traverse_data_t {
  int myRank;
  EW* ew;
  /* bool cartCoordSet; */
  string inFileName;
  string outFileName;
  int writeEvery;
  int downSample;
  TimeSeries::receiverMode mode;
  int event;
  vector< vector<TimeSeries*> > *GlobalTimeSeries;
  float_sw4 m_global_xmax;
  float_sw4 m_global_ymax;
  bool is_obs;
  bool winlset;
  bool winrset;
  float_sw4 winl;
  float_sw4 winr;
  float_sw4 win1;
  float_sw4 win2;
  float_sw4 win3;
  float_sw4 win4;
  bool usex;
  bool usey;
  bool usez;
  float_sw4 t0;
  bool scalefactor_set; 
  float_sw4 scalefactor;
} traverse_data_t;

struct traverse_data2_t {
  vector<string> *staname;
  vector<double> *x;
  vector<double> *y;
  vector<double> *z;
  vector<int>    *is_nsew;
  int            *n;
} traverse_data2_t;

struct srf_meta_t {
    float elon;
    float elat;
    int nstk;
    int ndip;
    float len;
    float wid;
    float stk;
    float dip;
    float dtop;
    float shyp;
    float dhyp;
} srf_meta_t;

struct srf_data_t {
    float lon;
    float lat;
    float dep;
    float stk;
    float dip;
    float area;
    float tinit;
    float dt;
    float vs;
    float den;
    float rake;
    float slip1;
    int   nt1;
    float slip2;
    int   nt2;
    float slip3;
    int   nt3;
} srf_data_t;

static herr_t traverse_func (hid_t loc_id, const char *grp_name, const H5L_info_t *info, void *operator_data)
{
  hid_t grp = -1, dset = -1, attr = -1;
  herr_t status;
#if H5_VERSION_GE(1,12,0)
  H5O_info1_t infobuf;
#else
  H5O_info_t infobuf;
#endif
  EW *a_ew;
  double data[4];
  double lon, lat, depth, x, y, z;
  bool geoCoordSet = true, topodepth = true, nsew = true;
  int isnsew, usezvalue, ret;
  bool foundwins=false;

  ASSERT(operator_data != NULL);

  struct traverse_data_t *op_data = (struct traverse_data_t *)operator_data;
  a_ew = op_data->ew;
  ASSERT(a_ew != NULL);

#if H5_VERSION_GE(1,12,0)
  status = H5Oget_info_by_name1(loc_id, grp_name, &infobuf, H5P_DEFAULT);
#else
  status = H5Oget_info_by_name(loc_id, grp_name, &infobuf, H5P_DEFAULT);
#endif
  if (infobuf.type == H5O_TYPE_GROUP) {
    /* if (op_data->myRank == 0) */
    /*   printf ("Group: [%s] \n", grp_name); */

    // read x,y,z or ns,ew,up
    grp = H5Gopen(loc_id, grp_name, H5P_DEFAULT);
    if (grp < 0) {
      fprintf(stderr, "Error opening group [%s]\n", grp_name);
      return -1;
    }

    if (H5Lexists(grp, "ISNSEW", H5P_DEFAULT) > 0) {
      attr = H5Dopen(grp, "ISNSEW", H5P_DEFAULT);
      ASSERT(attr > 0);
      ret = H5Dread(attr, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &isnsew);
      ASSERT(ret >= 0);
      H5Dclose(attr);
      if (isnsew == 0)
        nsew = false;
    }

    if (H5Lexists(grp, "WindowL", H5P_DEFAULT) > 0) {
      op_data->winlset = true;
      attr = H5Dopen(grp, "WindowL", H5P_DEFAULT);
      ASSERT(attr > 0);
      ret = H5Dread(attr, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &op_data->winl);
      ASSERT(ret >= 0);
      H5Dclose(attr);
    }

    if (H5Lexists(grp, "WindowR", H5P_DEFAULT) > 0) {
      op_data->winrset = true;
      attr = H5Dopen(grp, "WindowR", H5P_DEFAULT);
      ASSERT(attr > 0);
      ret = H5Dread(attr, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &op_data->winr);
      ASSERT(ret >= 0);
      H5Dclose(attr);
    }

    if (H5Lexists(grp, "USEZVALUE", H5P_DEFAULT) > 0) {
      attr = H5Dopen(grp, "USEZVALUE", H5P_DEFAULT);
      ASSERT(attr > 0);
      ret = H5Dread(attr, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &usezvalue);
      ASSERT(ret >= 0);
      H5Dclose(attr);
      if (usezvalue != 0)
        topodepth = false;
    }

    if (H5Lexists(grp, "WINDOWS", H5P_DEFAULT) > 0) {
       dset = H5Dopen(grp, "WINDOWS", H5P_DEFAULT);
       if (dset < 0)
          fprintf(stderr, "Error reading from rechdf5 station %s WINDOWS open failed!\n", grp_name);
       ret = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
       if( ret >= 0 )
          H5Dclose(dset);
       op_data->win1=data[0];
       op_data->win2=data[1];
       op_data->win3=data[2];
       op_data->win4=data[3];
       foundwins = true;
    }

    if (H5Lexists(grp, "STX,STY,STZ", H5P_DEFAULT) > 0) {
      // X, Y, Z
      dset = H5Dopen(grp, "STX,STY,STZ", H5P_DEFAULT);
      if (dset < 0)
        fprintf(stderr, "Error reading from rechdf5 station %s, STX,STY,STZ open failed!\n", grp_name);
      ASSERT(dset > 0);
      ret = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
      ASSERT(ret >= 0);
      H5Dclose(dset);
      x = data[0];
      y = data[1];
      z = data[2];
      geoCoordSet = false;
    }
    else if (H5Lexists(grp, "STLA,STLO,STDP", H5P_DEFAULT) > 0) {
      // STLA,STLO,STDP
      dset = H5Dopen(grp, "STLA,STLO,STDP", H5P_DEFAULT);
      if (dset < 0)
        fprintf(stderr, "Error reading from rechdf5 station %s, STLA,STLO,STDP open failed!\n", grp_name);
      ASSERT(dset > 0);
      ret = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
      ASSERT(ret >= 0);
      H5Dclose(dset);
      lat = data[0];
      lon = data[1];
      z = data[2];
    }
    else {
      // Not a station group, ignore
      H5Gclose(grp);
      return 0;
    }

    depth = z;
    if (geoCoordSet)
      a_ew->computeCartesianCoord(x, y, lon, lat);

    bool inCurvilinear=false;
    // we are in or above the curvilinear grid 
    if ( a_ew->topographyExists() && z < a_ew->m_zmin[a_ew->mNumberOfCartesianGrids-1])
      inCurvilinear = true;

    // check if (x,y,z) is not in the global bounding box
    if ( !( (inCurvilinear || z >= 0) && x>=0 && x<=op_data->m_global_xmax && y>=0 && y<=op_data->m_global_ymax)) {
      // The location of this station was outside the domain, so don't include it in the global list
      if (op_data->myRank == 0 && a_ew->getVerbosity() > 0) {
        stringstream receivererr;
    
        receivererr << endl 
  		  << "***************************************************" << endl
  		  << " WARNING:  RECEIVER positioned outside grid!" << endl;
        receivererr << " No RECEIVER file will be generated for file = " << op_data->outFileName<< endl;
        if (geoCoordSet) {
  	  receivererr << " @ lon=" << lon << " lat=" << lat << " depth=" << depth << endl << endl;
        }
        else {
    	  receivererr << " @ x=" << x << " y=" << y << " z=" << z << endl << endl;
        }
        
        receivererr << "***************************************************" << endl;
        cerr << receivererr.str();
        cerr.flush();
      }
    }
    else
    {
      /* if (op_data->myRank == 0) */
      /*   cout << "x=" << x << ", y=" << y << ", z=" << z << ", writeEvery=" << op_data->writeEvery << endl; */

      TimeSeries *ts_ptr = new TimeSeries(a_ew, grp_name, grp_name, op_data->mode, false, false, true, op_data->outFileName, x, y, z, 
  					topodepth, op_data->writeEvery, op_data->downSample, !nsew, op_data->event );

      if((*op_data->GlobalTimeSeries)[op_data->event].size() == 0) {
        ts_ptr->allocFid();
        ts_ptr->setTS0Ptr(ts_ptr);
      }
      else {
        ts_ptr->setFidPtr((*op_data->GlobalTimeSeries)[op_data->event][0]->getFidPtr());
        ts_ptr->setTS0Ptr((*op_data->GlobalTimeSeries)[op_data->event][0]);
      }
     
      if (ts_ptr->myPoint()) {
        /* cout << "Rank " << op_data->myRank << "has this point x=" << x << " y=" << y << " z=" << z << endl; */

        // Only for observation data
        if (op_data->is_obs) {
          // Read data
          bool ignore_utc = false;
          ts_ptr->readSACHDF5(op_data->ew, op_data->inFileName, ignore_utc);

          // Set reference UTC to simulation UTC, for easier plotting.
          ts_ptr->set_utc_to_simulation_utc();
      
          // Set window, in simulation time
          if( op_data->winlset || op_data->winrset )
          {
             if( op_data->winlset && !op_data->winrset )
                op_data->winr = 1e38;
             if( !op_data->winlset && op_data->winrset )
                op_data->winl = -1;
             ts_ptr->set_window( op_data->winl, op_data->winr );
          }

          if( foundwins )
             ts_ptr->set_window( op_data->win1, op_data->win2, op_data->win3, op_data->win4 );

          // Exclude some components
          if( !op_data->usex || !op_data->usey || !op_data->usez )
             ts_ptr->exclude_component( op_data->usex, op_data->usey, op_data->usez );
      
          // Add extra shift from command line, use with care.
          if( op_data->t0 != 0 )
             ts_ptr->add_shift( op_data->t0 );
          //DBG
          //          ts_ptr->set_shift(0.0);
      
          // Set scale factor if given
          if( op_data->scalefactor_set )
             ts_ptr->set_scalefactor( op_data->scalefactor );
          }
      }
  
      // include the receiver in the global list
      (*op_data->GlobalTimeSeries)[op_data->event].push_back(ts_ptr);
    }

    H5Gclose(grp);
  }

  return 0;
}

// Callback for the cheap "name only" scan over the station file's root group.
// It records every link name without opening the target object, so it does no
// per-object metadata reads - the expensive per-station dataset reads are
// deferred to the parallel phase in readStationHDF5 below. Runs on rank 0 only.
static herr_t collect_names_cb(hid_t loc_id, const char *grp_name, const H5L_info_t *info, void *operator_data)
{
  vector<string> *names = (vector<string> *)operator_data;
  names->push_back(grp_name);
  return 0;
}

// One station's data as read from the file.
struct sta_rec_t { double x, y, z; int nsew; int topodepth; };

// Read a single station group by name and fill 'out'. Returns true if 'name' is
// a valid station group whose location is inside the domain bounding box
// [0,gxmax] x [0,gymax] (matching the original per-rank traversal's accept
// test). Prints the "outside grid" warning for rejected stations when verbose.
// Pure per-rank work: no MPI, so it can run concurrently on many reader ranks,
// each over a disjoint set of names.
static bool read_one_station(hid_t file_id, const string &name, EW *ew,
                             double gxmax, double gymax, sta_rec_t &out)
{
  double data[4];
  double lon = 0, lat = 0, depth, x = 0, y = 0, z = 0;
  bool geoCoordSet = true, topodepth = true, nsew = true;
  int isnsew, usezvalue, ret;

  hid_t grp;
  H5E_BEGIN_TRY { grp = H5Gopen(file_id, name.c_str(), H5P_DEFAULT); } H5E_END_TRY;
  if (grp < 0)
    return false;   // link is not an openable group

  if (H5Lexists(grp, "ISNSEW", H5P_DEFAULT) > 0) {
    hid_t attr = H5Dopen(grp, "ISNSEW", H5P_DEFAULT);
    ret = H5Dread(attr, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &isnsew);
    H5Dclose(attr);
    if (ret >= 0 && isnsew == 0)
      nsew = false;
  }

  if (H5Lexists(grp, "USEZVALUE", H5P_DEFAULT) > 0) {
    hid_t attr = H5Dopen(grp, "USEZVALUE", H5P_DEFAULT);
    ret = H5Dread(attr, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &usezvalue);
    H5Dclose(attr);
    if (ret >= 0 && usezvalue != 0)
      topodepth = false;
  }

  if (H5Lexists(grp, "STX,STY,STZ", H5P_DEFAULT) > 0) {
    hid_t dset = H5Dopen(grp, "STX,STY,STZ", H5P_DEFAULT);
    if (dset < 0) {
      fprintf(stderr, "Error reading from rechdf5 station %s, STX,STY,STZ open failed!\n", name.c_str());
      H5Gclose(grp);
      return false;
    }
    ret = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(dset);
    x = data[0];
    y = data[1];
    z = data[2];
    geoCoordSet = false;
  }
  else if (H5Lexists(grp, "STLA,STLO,STDP", H5P_DEFAULT) > 0) {
    hid_t dset = H5Dopen(grp, "STLA,STLO,STDP", H5P_DEFAULT);
    if (dset < 0) {
      fprintf(stderr, "Error reading from rechdf5 station %s, STLA,STLO,STDP open failed!\n", name.c_str());
      H5Gclose(grp);
      return false;
    }
    ret = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(dset);
    lat = data[0];
    lon = data[1];
    z = data[2];
  }
  else {
    // Not a station group, ignore
    H5Gclose(grp);
    return false;
  }
  H5Gclose(grp);

  depth = z;
  if (geoCoordSet)
    ew->computeCartesianCoord(x, y, lon, lat);

  bool inCurvilinear = false;
  if (ew->topographyExists() && z < ew->m_zmin[ew->mNumberOfCartesianGrids-1])
    inCurvilinear = true;

  if (!((inCurvilinear || z >= 0) && x >= 0 && x <= gxmax && y >= 0 && y <= gymax)) {
    // Outside the domain: don't include it in the global list
    if (ew->getVerbosity() > 0) {
      stringstream receivererr;
      receivererr << endl
                  << "***************************************************" << endl
                  << " WARNING:  RECEIVER positioned outside grid!" << endl
                  << " No RECEIVER file will be generated for station = " << name << endl;
      if (geoCoordSet)
        receivererr << " @ lon=" << lon << " lat=" << lat << " depth=" << depth << endl << endl;
      else
        receivererr << " @ x=" << x << " y=" << y << " z=" << z << endl << endl;
      receivererr << "***************************************************" << endl;
      cerr << receivererr.str();
      cerr.flush();
    }
    return false;
  }

  out.x = x;
  out.y = y;
  out.z = z;
  out.nsew = nsew ? 1 : 0;
  out.topodepth = topodepth ? 1 : 0;
  return true;
}

void readStationHDF5(EW* ew, string inFileName, string outFileName, int writeEvery, int downSample, TimeSeries::receiverMode mode, int event, vector< vector<TimeSeries*> > *GlobalTimeSeries, float_sw4 m_global_xmax, float_sw4 m_global_ymax, bool is_obs, bool winlset, bool winrset, float_sw4 winl, float_sw4 winr, bool usex, bool usey, bool usez, float_sw4 t0, bool scalefactor_set, float_sw4 scalefactor)
{
  if (!is_obs) {
    // Bulk receiver path (rechdf5/receiver command). Built to scale to 1e5+
    // stations. The original code had every rank independently H5Literate the
    // whole file (~10 HDF5 metadata reads per station), i.e. O(stations x ranks)
    // uncoordinated filesystem traffic, plus 2 collective MPI_Allreduce per
    // station inside the TimeSeries constructor, i.e. O(stations) serialized
    // collectives. Both dominate startup at scale. This replacement is:
    //
    //   1. rank 0 enumerates station group names only (cheap link iteration,
    //      no per-object reads) and broadcasts them.
    //   2. a spread-out subset of "reader" ranks each open the file once and
    //      read a disjoint block of stations in parallel; results are
    //      Allgatherv'd to all ranks. Total dataset reads = O(stations).
    //   3. topography elevation and the "exactly one owner" safety check are
    //      each done as a single batched Allreduce over all stations, instead
    //      of per-station collectives (see TimeSeries ctor deferCollectives).
    MPI_Comm comm = ew->m_1d_communicator;
    int rank = ew->getRank();
    int nprocs = 1;
    MPI_Comm_size(comm, &nprocs);
    double t0all = MPI_Wtime();

    // ---- 1. rank 0 enumerates station names, broadcast to all ranks ----
    vector<string> names;
    if (rank == 0) {
      if (inFileName == outFileName)
        printf("Warning: Same station input file and output file name [%s]\n", inFileName.c_str());
      printf("readStationHDF5: enumerating stations in %s ...\n", inFileName.c_str());
      fflush(stdout);

      hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
      hid_t fid = H5Fopen(inFileName.c_str(), H5F_ACC_RDONLY, fapl);
      if (fid < 0)
        printf("%s Error opening file [%s]\n", __func__, inFileName.c_str());
      else {
        H5Literate(fid, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, collect_names_cb, &names);
        H5Fclose(fid);
      }
      H5Pclose(fapl);

      printf("readStationHDF5: found %d candidate station groups in %.2f s\n",
             (int)names.size(), MPI_Wtime() - t0all);
      fflush(stdout);
    }

    int ncand = (int)names.size();
    MPI_Bcast(&ncand, 1, MPI_INT, 0, comm);

    // Broadcast the names ('\0'-separated buffer).
    string namebuf;
    if (rank == 0)
      for (int s = 0; s < ncand; s++) { namebuf += names[s]; namebuf += '\0'; }
    int namebuflen = (int)namebuf.size();
    MPI_Bcast(&namebuflen, 1, MPI_INT, 0, comm);
    vector<char> namebytes(namebuflen > 0 ? namebuflen : 1);
    if (rank == 0 && namebuflen > 0)
      memcpy(namebytes.data(), namebuf.data(), namebuflen);
    if (namebuflen > 0)
      MPI_Bcast(namebytes.data(), namebuflen, MPI_CHAR, 0, comm);
    if (rank != 0) {
      names.resize(ncand);
      int pos = 0;
      for (int s = 0; s < ncand; s++) { names[s] = string(&namebytes[pos]); pos += (int)names[s].size() + 1; }
    }

    if (ncand == 0)
      return;

    // ---- 2. parallel read of station data by a spread-out reader subset ----
    // Cap the number of concurrent readers so we don't overload the filesystem
    // metadata server; spread them evenly across the job (every 'stride'-th rank)
    // so the readers land on as many different nodes as possible.
    int target_readers = nprocs < 128 ? nprocs : 128;
    if (target_readers > ncand) target_readers = ncand;
    if (target_readers < 1) target_readers = 1;
    int stride = nprocs / target_readers;              // >= 1
    if (stride < 1) stride = 1;
    bool am_reader = (rank % stride == 0) && (rank / stride < target_readers);
    int reader_id  = am_reader ? rank / stride : -1;
    int nreaders   = target_readers;                   // reader_id in [0, nreaders)
    int block = (ncand + nreaders - 1) / nreaders;     // names per reader

    vector<double> lcoord;    // accepted x,y,z triples read by this rank
    vector<int>    lflag;     // accepted nsew,topodepth pairs read by this rank
    string         lnames;    // accepted names, '\0'-separated
    double t0read = MPI_Wtime();
    if (am_reader) {
      int begin = reader_id * block;
      int end   = begin + block;
      if (end > ncand) end = ncand;
      hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
      hid_t fid  = H5Fopen(inFileName.c_str(), H5F_ACC_RDONLY, fapl);
      if (fid < 0)
        printf("%s reader rank %d: error opening file [%s]\n", __func__, rank, inFileName.c_str());
      else {
        for (int s = begin; s < end; s++) {
          sta_rec_t rec;
          if (read_one_station(fid, names[s], ew, m_global_xmax, m_global_ymax, rec)) {
            lcoord.push_back(rec.x); lcoord.push_back(rec.y); lcoord.push_back(rec.z);
            lflag.push_back(rec.nsew); lflag.push_back(rec.topodepth);
            lnames += names[s]; lnames += '\0';
          }
        }
        H5Fclose(fid);
      }
      H5Pclose(fapl);
    }

    // ---- Allgatherv the accepted stations to all ranks ----
    // Every rank ends up with an identical global list (reader-block order),
    // so the batched Allreduces below index the same station on every rank.
    int local_n = (int)lflag.size() / 2;
    vector<int> counts(nprocs), displs(nprocs);
    MPI_Allgather(&local_n, 1, MPI_INT, counts.data(), 1, MPI_INT, comm);
    int total = 0;
    for (int p = 0; p < nprocs; p++) { displs[p] = total; total += counts[p]; }

    if (rank == 0) {
      printf("readStationHDF5: %d reader ranks read %d in-domain stations in %.2f s\n",
             nreaders, total, MPI_Wtime() - t0read);
      fflush(stdout);
    }
    if (total == 0)
      return;

    // coordinates: 3 doubles per station
    vector<int> ccounts(nprocs), cdispls(nprocs);
    for (int p = 0; p < nprocs; p++) { ccounts[p] = 3*counts[p]; cdispls[p] = 3*displs[p]; }
    vector<double> gcoord(3*total);
    MPI_Allgatherv(lcoord.data(), 3*local_n, MPI_DOUBLE,
                   gcoord.data(), ccounts.data(), cdispls.data(), MPI_DOUBLE, comm);
    // flags: 2 ints per station
    vector<int> fcounts(nprocs), fdispls(nprocs);
    for (int p = 0; p < nprocs; p++) { fcounts[p] = 2*counts[p]; fdispls[p] = 2*displs[p]; }
    vector<int> gflag(2*total);
    MPI_Allgatherv(lflag.data(), 2*local_n, MPI_INT,
                   gflag.data(), fcounts.data(), fdispls.data(), MPI_INT, comm);
    // names: '\0'-separated bytes
    int lnbytes = (int)lnames.size();
    vector<int> ncounts(nprocs), ndispls(nprocs);
    MPI_Allgather(&lnbytes, 1, MPI_INT, ncounts.data(), 1, MPI_INT, comm);
    int totbytes = 0;
    for (int p = 0; p < nprocs; p++) { ndispls[p] = totbytes; totbytes += ncounts[p]; }
    vector<char> gnames(totbytes > 0 ? totbytes : 1);
    MPI_Allgatherv(lnames.data(), lnbytes, MPI_CHAR,
                   gnames.data(), ncounts.data(), ndispls.data(), MPI_CHAR, comm);
    vector<string> snames(total);
    { int pos = 0; for (int s = 0; s < total; s++) { snames[s] = string(&gnames[pos]); pos += (int)snames[s].size() + 1; } }

    // ---- 3a. batched topography elevation: one Allreduce over all stations ----
    // Each rank interpolates topography locally for every station (returning a
    // sentinel where the (x,y) column isn't in its subdomain); one MPI_MAX
    // reduction then yields the true elevation, replacing the per-station
    // Allreduce the TimeSeries ctor would otherwise do.
    vector<double> zTopo(total, 0.0);
    if (ew->topographyExists()) {
      vector<float_sw4> zloc(total), zglob(total);
      for (int s = 0; s < total; s++) {
        float_sw4 zt;
        if (!ew->m_gridGenerator->interpolate_topography(ew, gcoord[3*s], gcoord[3*s+1], zt, ew->mTopoGridExt))
          zt = -1e38;
        zloc[s] = zt;
      }
      MPI_Allreduce(zloc.data(), zglob.data(), total, ew->m_mpifloat, MPI_MAX, comm);
      for (int s = 0; s < total; s++) zTopo[s] = zglob[s];
    }

    // ---- 3b. construct TimeSeries objects (no per-station collectives) ----
    double t0con = MPI_Wtime();
    int base = (int)(*GlobalTimeSeries)[event].size();
    for (int s = 0; s < total; s++) {
      TimeSeries *ts_ptr = new TimeSeries(ew, snames[s], snames[s], mode, false, false, true,
                                          outFileName, gcoord[3*s], gcoord[3*s+1], gcoord[3*s+2],
                                          gflag[2*s+1] != 0, writeEvery, downSample,
                                          gflag[2*s] == 0, event,
                                          /*deferCollectives=*/true, (float_sw4)zTopo[s]);
      if ((*GlobalTimeSeries)[event].size() == 0) {
        ts_ptr->allocFid();
        ts_ptr->setTS0Ptr(ts_ptr);
      }
      else {
        ts_ptr->setFidPtr((*GlobalTimeSeries)[event][0]->getFidPtr());
        ts_ptr->setTS0Ptr((*GlobalTimeSeries)[event][0]);
      }
      (*GlobalTimeSeries)[event].push_back(ts_ptr);
    }

    // ---- 3c. batched "exactly one owner" safety check: one Allreduce ----
    {
      vector<int> owned(total), ownsum(total);
      for (int s = 0; s < total; s++)
        owned[s] = (*GlobalTimeSeries)[event][base + s]->myPoint() ? 1 : 0;
      MPI_Allreduce(owned.data(), ownsum.data(), total, MPI_INT, MPI_SUM, comm);
      if (rank == 0) {
        int nbad = 0;
        for (int s = 0; s < total; s++) if (ownsum[s] != 1) nbad++;
        if (nbad > 0)
          printf("readStationHDF5: WARNING, %d of %d stations are not owned by exactly one rank\n",
                 nbad, total);
      }
    }

    if (rank == 0) {
      printf("readStationHDF5: constructed %d TimeSeries in %.2f s (total station setup %.2f s)\n",
             total, MPI_Wtime() - t0con, MPI_Wtime() - t0all);
      fflush(stdout);
    }
    return;
  }

  // is_obs == true (observation command): each owning rank must read its
  // own SAC/HDF5 waveform data out of the same file via readSACHDF5(), so
  // every rank still needs its own file handle - left as a per-rank
  // H5Literate traversal.
  hid_t fid, fapl;

  struct traverse_data_t tData;
  tData.myRank = ew->getRank();
  tData.ew = ew;
  tData.inFileName  = inFileName;
  tData.outFileName = outFileName;
  tData.writeEvery  = writeEvery;
  tData.downSample  = downSample;
  tData.mode = mode;
  tData.event = event;
  tData.m_global_xmax = m_global_xmax;
  tData.m_global_ymax = m_global_ymax;
  tData.GlobalTimeSeries = GlobalTimeSeries;
  tData.is_obs = is_obs;
  tData.winlset = winlset;
  tData.winrset = winrset;
  tData.winl = winl;
  tData.winr = winr;
  tData.usex = usex;
  tData.usey = usey;
  tData.usez = usez;
  tData.t0 = t0;
  tData.scalefactor_set = scalefactor_set;
  tData.scalefactor = scalefactor;

  fapl = H5Pcreate(H5P_FILE_ACCESS);
  /* H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL); */

  fid = H5Fopen(inFileName.c_str(),  H5F_ACC_RDONLY, fapl);
  if (fid < 0) {
    printf("%s Error opening file [%s]\n", __func__, inFileName.c_str());
    return;
  }

  if (inFileName == outFileName) {
    if (tData.myRank == 0) {
      printf("Warning: Same station input file and output file name [%s]\n", inFileName.c_str());
    }
  }

  H5Literate (fid, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, traverse_func, &tData);

  H5Pclose(fapl);
  H5Fclose(fid);

  return;
}

static herr_t traverse_func2 (hid_t loc_id, const char *grp_name, const H5L_info_t *info, void *operator_data)
{
  hid_t grp = -1, dset = -1, attr = -1;
  herr_t status;
#if H5_VERSION_GE(1,12,0)
  H5O_info1_t infobuf;
#else
  H5O_info_t infobuf;
#endif
  float data[3];
  int isnsew, ret;

  ASSERT(operator_data != NULL);

  struct traverse_data2_t *op_data = (struct traverse_data2_t *)operator_data;


#if H5_VERSION_GE(1,12,0)
  status = H5Oget_info_by_name1(loc_id, grp_name, &infobuf, H5P_DEFAULT);
#else
  status = H5Oget_info_by_name(loc_id, grp_name, &infobuf, H5P_DEFAULT);
#endif
  if (infobuf.type == H5O_TYPE_GROUP) {
    /* if (op_data->myRank == 0) */
    /*   printf ("Group: [%s] \n", grp_name); */

    // read x,y,z or ns,ew,up
    grp = H5Gopen(loc_id, grp_name, H5P_DEFAULT);
    if (grp < 0) {
      printf("Error opening group [%s]\n", grp_name);
      return -1;
    }

    if (H5Lexists(grp, "ISNSEW", H5P_DEFAULT) > 0) {
      attr = H5Dopen(grp, "ISNSEW", H5P_DEFAULT);
      ASSERT(attr > 0);
      ret = H5Dread(attr, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &isnsew);
      ASSERT(ret >= 0);
      H5Dclose(attr);
    }

    if (H5Lexists(grp, "STX,STY,STZ", H5P_DEFAULT) > 0) {
      // X, Y, Z
      dset = H5Dopen(grp, "STX,STY,STZ", H5P_DEFAULT);
      if (dset < 0)
        fprintf(stderr, "Error reading from rechdf5 station %s, STX,STY,STZ open failed!\n", grp_name);
      ASSERT(dset > 0);
      ret = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
      ASSERT(ret >= 0);
      H5Dclose(dset);
    }
    else if (H5Lexists(grp, "STLA,STLO,STDP", H5P_DEFAULT) > 0) {
      // STLA,STLO,STDP
      dset = H5Dopen(grp, "STLA,STLO,STDP", H5P_DEFAULT);
      if (dset < 0)
        fprintf(stderr, "Error reading from rechdf5 station %s, STLA,STLO,STDP open failed!\n", grp_name);
      ASSERT(dset > 0);
      ret = H5Dread(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
      ASSERT(ret >= 0);
      H5Dclose(dset);
    }
    else {
      // Not a station group, ignore
      H5Gclose(grp);
      return 0;
    }

    string staname = grp_name;
    (*op_data->staname).push_back(staname);
    (*op_data->x).push_back(data[0]);
    (*op_data->y).push_back(data[1]);
    (*op_data->z).push_back(data[2]);
    (*op_data->is_nsew).push_back(isnsew);
    (*op_data->n)++;

    H5Gclose(grp);
  }

  return 0;
}

void readStationInfoHDF5(string inFileName, vector<string> *staname, vector<double> *x, vector<double> *y, vector<double> *z, vector<int> *is_nsew, int *n)
{
  hid_t fid, fapl;

  struct traverse_data2_t tData;
  tData.staname = staname;
  tData.x = x;
  tData.y = y;
  tData.z = z;
  tData.is_nsew = is_nsew;
  tData.n = n;

  fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fapl, MPI_COMM_SELF, MPI_INFO_NULL);

  fid = H5Fopen(inFileName.c_str(),  H5F_ACC_RDONLY, fapl);
  if (fid < 0) {
    printf("%s Error opening file [%s]\n", __func__, inFileName.c_str());
    return;
  }

  H5Literate (fid, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, traverse_func2, &tData);

  H5Pclose(fapl);
  H5Fclose(fid);

  return;
}

void readRuptureHDF5(char *fname, vector<vector<Source*> > & a_GlobalUniqueSources, EW *ew, int event, float_sw4 m_global_xmax, float_sw4 m_global_ymax, float_sw4 m_global_zmax, float_sw4 mGeoAz, float_sw4 xmin, float_sw4 ymin, float_sw4 zmin, int mVerbose, int nreader)
{
  bool is_debug = true;
  int world_rank, world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  if (nreader <= 0) 
      nreader = 1;
  if (nreader > world_size)
      nreader = world_size;

  int read_color = world_rank % (world_size / nreader) == 0 ? 0 : 1;
  int node_color = world_rank / (world_size / nreader);
  int read_rank, read_size;
  MPI_Comm read_comm, node_comm;
  MPI_Comm_split(MPI_COMM_WORLD, read_color, world_rank, &read_comm);
  MPI_Comm_split(MPI_COMM_WORLD, node_color, world_rank, &node_comm);
  MPI_Comm_rank(MPI_COMM_WORLD, &read_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &read_size);

  double stime, etime;
  hid_t fid, grp, attr, ctype, dtype, dset, dspace, aspace, fapl;

  ctype = H5Tcreate(H5T_COMPOUND, 9 * sizeof(float) + 2 * sizeof(int));
  H5Tinsert(ctype, "ELON",                  0, H5T_NATIVE_FLOAT);
  H5Tinsert(ctype, "ELAT",  1 * sizeof(float), H5T_NATIVE_FLOAT);
  H5Tinsert(ctype, "NSTK",  2 * sizeof(float), H5T_NATIVE_INT);
  H5Tinsert(ctype, "NDIP",  2 * sizeof(float) + 1 * sizeof(int), H5T_NATIVE_INT);
  H5Tinsert(ctype, "LEN",   2 * sizeof(float) + 2 * sizeof(int), H5T_NATIVE_FLOAT);
  H5Tinsert(ctype, "WID",   3 * sizeof(float) + 2 * sizeof(int), H5T_NATIVE_FLOAT);
  H5Tinsert(ctype, "STK",   4 * sizeof(float) + 2 * sizeof(int), H5T_NATIVE_FLOAT);
  H5Tinsert(ctype, "DIP",   5 * sizeof(float) + 2 * sizeof(int), H5T_NATIVE_FLOAT);
  H5Tinsert(ctype, "DTOP",  6 * sizeof(float) + 2 * sizeof(int), H5T_NATIVE_FLOAT);
  H5Tinsert(ctype, "SHYP",  7 * sizeof(float) + 2 * sizeof(int), H5T_NATIVE_FLOAT);
  H5Tinsert(ctype, "DHYP",  8 * sizeof(float) + 2 * sizeof(int), H5T_NATIVE_FLOAT);

  dtype = H5Tcreate(H5T_COMPOUND, 14 * sizeof(float) + 3 * sizeof(int));
  H5Tinsert(dtype, "LON", 0, H5T_NATIVE_FLOAT);
  H5Tinsert(dtype, "LAT", 1 * sizeof(float), H5T_NATIVE_FLOAT);
  H5Tinsert(dtype, "DEP", 2 * sizeof(float), H5T_NATIVE_FLOAT);
  H5Tinsert(dtype, "STK", 3 * sizeof(float), H5T_NATIVE_FLOAT);
  H5Tinsert(dtype, "DIP", 4 * sizeof(float), H5T_NATIVE_FLOAT);
  H5Tinsert(dtype, "AREA", 5 * sizeof(float), H5T_NATIVE_FLOAT);
  H5Tinsert(dtype, "TINIT", 6 * sizeof(float), H5T_NATIVE_FLOAT);
  H5Tinsert(dtype, "DT", 7 * sizeof(float), H5T_NATIVE_FLOAT);
  H5Tinsert(dtype, "VS", 8 * sizeof(float), H5T_NATIVE_FLOAT);
  H5Tinsert(dtype, "DEN", 9 * sizeof(float), H5T_NATIVE_FLOAT);
  H5Tinsert(dtype, "RAKE", 10 * sizeof(float), H5T_NATIVE_FLOAT);
  H5Tinsert(dtype, "SLIP1", 11 * sizeof(float), H5T_NATIVE_FLOAT);
  H5Tinsert(dtype, "NT1", 12 * sizeof(float), H5T_NATIVE_INT);
  H5Tinsert(dtype, "SLIP2", 12 * sizeof(float) + 1 * sizeof(int), H5T_NATIVE_FLOAT);
  H5Tinsert(dtype, "NT2", 13 * sizeof(float) + 1 * sizeof(int), H5T_NATIVE_INT);
  H5Tinsert(dtype, "SLIP3", 13 * sizeof(float) + 2 * sizeof(int), H5T_NATIVE_FLOAT);
  H5Tinsert(dtype, "NT3", 14 * sizeof(float) + 2 * sizeof(int), H5T_NATIVE_INT);

  struct srf_meta_t *srf_metadata;
  struct srf_data_t *point_data;
  float *sr_data;

  int npts = 0, nseg = 0, nsr1 = 0;
  hsize_t dims;
  double rVersion;
  int nSources=0, nu1=0, nu2=0, nu3=0, nskip_zero_slip=0;

  stime = MPI_Wtime();
  // Only rank 0 reads data, then broadcast to all other processes
  if (read_color == 0) {

    fapl = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(fapl, read_comm, MPI_INFO_NULL);
    fid = H5Fopen(fname, H5F_ACC_RDONLY, fapl);
    if (fid <= 0) 
      cout << "Rupture HDF5 file " << fname << " not found" << endl;
    
    if (world_rank == 0) 
      printf("Opened rupture file '%s'\n", fname);

    attr = H5Aopen(fid, "VERSION", H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_DOUBLE, &rVersion);
    H5Aclose(attr);
    if (world_rank == 0) 
      printf("Version = %.1f\n", rVersion);

    // read each header block
    attr = H5Aopen(fid, "PLANE", H5P_DEFAULT);
    aspace = H5Aget_space(attr);
    H5Sget_simple_extent_dims(aspace, &dims, NULL);
    nseg = (int)dims;
    srf_metadata = (struct srf_meta_t *)malloc(nseg * sizeof(struct srf_meta_t));
    H5Sclose(aspace);
    if (world_rank == 0) 
      printf("Number of segments in header block: %i\n", nseg);
    H5Aread(attr, ctype, srf_metadata);
    H5Aclose(attr);

    if (world_rank == 0) {
      for (int seg=0; seg<nseg; seg++) {
        printf("Seg #%i: elon=%g, elat=%g, nstk=%i, ndip=%i, len=%g, wid=%g\n", 
                seg+1, srf_metadata[seg].elon, srf_metadata[seg].elat, srf_metadata[seg].nstk, srf_metadata[seg].ndip, srf_metadata[seg].len, srf_metadata[seg].wid);
        printf("        stk=%g, dip=%g, dtop=%g, shyp=%g, dhyp=%g\n", 
                srf_metadata[seg].stk, srf_metadata[seg].dip, srf_metadata[seg].dtop, srf_metadata[seg].shyp, srf_metadata[seg].dhyp);
      }
    }

    free(srf_metadata);
    dset = H5Dopen(fid, "POINTS", H5P_DEFAULT);
    if (dset < 0) {
      printf("Error with Rupture HDF5 file, no POINTS dataset found!\n");
      npts = -1;
    }
    else {
      dspace = H5Dget_space(dset);
 
      H5Sget_simple_extent_dims(dspace, &dims, NULL);
      npts = (int)dims;
      if (world_rank == 0) 
        printf("Number of point sources in data block: %i\n", npts);

      point_data = (struct srf_data_t*)malloc(npts*sizeof(struct srf_data_t));
      H5Dread(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, point_data);
      H5Sclose(dspace);
      H5Dclose(dset);
    }

    dset = H5Dopen(fid, "SR1", H5P_DEFAULT);
    if (dset < 0) {
      printf("Error with Rupture HDF5 file, no SR1 dataset found!\n");
      nsr1 = -1;
    }
    else {
      dspace = H5Dget_space(dset);
 
      H5Sget_simple_extent_dims(dspace, &dims, NULL);
      nsr1 = (int)dims;

      sr_data = (float*)malloc(nsr1*sizeof(float));
      H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, sr_data);
      H5Sclose(dspace);
      H5Dclose(dset);
    }

    H5Fclose(fid);
  }// End read_color=0
  etime = MPI_Wtime();

  if (is_debug && world_rank == 0) 
      printf("Read SRF-HDF5 takes %.2f seconds\n", etime-stime);

  MPI_Bcast(&npts, 1, MPI_INT, 0, node_comm);
  /* MPI_Bcast(&npts, 1, MPI_INT, 0, MPI_COMM_WORLD); */
  if (npts == -1) 
    return;

  MPI_Bcast(&nsr1, 1, MPI_INT, 0, node_comm);
  /* MPI_Bcast(&nsr1, 1, MPI_INT, 0, MPI_COMM_WORLD); */
  if (nsr1 == -1) 
    return;

  if (read_color != 0) {
    point_data = (struct srf_data_t*)malloc(npts*sizeof(struct srf_data_t));
    sr_data = (float*)malloc(nsr1*sizeof(float));
  }

  MPI_Bcast(point_data, npts*sizeof(struct srf_data_t), MPI_CHAR, 0, node_comm);
  /* MPI_Bcast(point_data, npts*sizeof(struct srf_data_t), MPI_CHAR, 0, MPI_COMM_WORLD); */

  MPI_Bcast(sr_data, nsr1, MPI_FLOAT, 0, node_comm);
  /* MPI_Bcast(sr_data, nsr1, MPI_FLOAT, 0, MPI_COMM_WORLD); */
  stime = MPI_Wtime();

  MPI_Comm_free(&node_comm);
  MPI_Comm_free(&read_comm);

  if (is_debug && world_rank == 0) 
      printf("Bcast SRF-HDF5 takes %.2f seconds\n", stime-etime);

  Source* sourcePtr;
  timeDep tDep = iDiscrete;
  char formstring[100];
  strcpy(formstring, "Discrete");

  double x = 0.0, y = 0.0, z = 0.0;
  float_sw4 m0 = 1.0;
  float_sw4 t0=0.0, f0=1.0, freq=1.0;
  float_sw4 mxx=0.0, mxy=0.0, mxz=0.0, myy=0.0, myz=0.0, mzz=0.0;
  int isMomentType = -1;
  bool topodepth = true;
  // Discrete source time function
  float_sw4* par=NULL;
  int* ipar=NULL;
  int npar=0, nipar=0, ncyc=0, sr1pos=0;
  // read all point sources
  for (int pts=0; pts<npts; pts++) 
  {
    double lon, lat, dep, stk, dip, area, tinit, dt, rake, slip1, slip2, slip3;
    int nt1, nt2, nt3;

    lon = (double)point_data[pts].lon;
    lat = (double)point_data[pts].lat;
    dep = (double)point_data[pts].dep;
    stk = (double)point_data[pts].stk;
    dip = (double)point_data[pts].dip;
    area = (double)point_data[pts].area;
    tinit = (double)point_data[pts].tinit;
    dt = (double)point_data[pts].dt;
    rake = (double)point_data[pts].rake;
    slip1 = (double)point_data[pts].slip1;
    slip2 = (double)point_data[pts].slip2;
    slip3 = (double)point_data[pts].slip3;
    nt1 = (int)point_data[pts].nt1;
    nt2 = (int)point_data[pts].nt2;
    nt3 = (int)point_data[pts].nt3;

    // nothing to do if nt1=nt2=nt3=0
    if (nt1<=0 && nt2<=0 && nt3<=0) continue;

    if (world_rank == 0 && mVerbose >= 2)
    {
      printf("point #%i: lon=%g, lat=%g, dep=%g, stk=%g, dip=%g, area=%g, tinit=%g, dt=%g\n", 
             pts+1, lon, lat, dep, stk, dip, area, tinit, dt);
      printf("          rake=%g, slip1=%g, nt1=%i, slip2=%g, nt2=%i, slip3=%g, nt3=%i\n", 
             rake, slip1, nt1, slip2, nt2, slip3, nt3);
    }
    
  // read discrete time series for u1
    if (nt1>0)
    {
      nu1++;
      // note that the first data point is always zero, but the last is not
      // for this reason we always pad the time zeries with a '0' 
      // also note that we need at least 7 data points, i.e. nt1>=6
      int nt1dim = max(6,nt1);
      par = new float_sw4[nt1dim+2];
      par[0]  = tinit;
      t0      = tinit;
      freq    = 1/dt;
      ipar    = new int[1];
      ipar[0] = nt1dim+1; // add an extra point 

      for( int i=0 ; i < nt1 ; i++ ) 
        par[i+1] = sr_data[sr1pos++];
      
      // pad with 0
      if (nt1 < 6) {
        for (int j=nt1; j<6; j++)
          par[j+1]=0.;
      }
      
      // last 0
      par[nt1dim+1]= 0.0;
  
      // scale cm/s to m/s
      for (int i=1; i<=nt1dim+1; i++)
      {
        par[i] *= 1e-2;
      }
  
      // AP: Mar. 1, 2016: Additional scaling is needed to make the integral of the time function = 1
      float_sw4 slip_m=slip1*1e-2;
      float_sw4 slip_sum=0;
      for (int i=1; i<=nt1dim+1; i++)
      {
        slip_sum += par[i];
      }
      slip_sum *=dt;
  
      if (world_rank == 0 && mVerbose >= 3)
      {
         printf("INFO: SRF file: dt*sum(slip_vel)=%e [m], total slip (from header)=%e [m]\n", slip_sum, slip_m);
      }
      float_sw4 slip_sum_tol = 1e-12;
      bool skip_zero_slip_point = false;
      if( slip_sum > -slip_sum_tol && slip_sum < slip_sum_tol )
      {
        nskip_zero_slip++;
        skip_zero_slip_point = true;
        if( world_rank == 0 && nskip_zero_slip <= 10 )
        {
          printf("WARNING: skipping rupture point #%i because dt*sum(slip_vel)=%e [m], total slip (from header)=%e [m]\n",
                 pts+1, slip_sum, slip_m);
        }
      }
      // scale time series to sum to integrate to one
      if( !skip_zero_slip_point )
      {
        for (int i=1; i<=nt1dim+1; i++)
        {
           par[i] /= slip_sum;
        }
        if (world_rank == 0 && mVerbose >= 3)
        {
           slip_sum=0;
           for (int i=1; i<=nt1dim+1; i++)
           {
              slip_sum += par[i];
           }
           slip_sum *=dt;
           printf("INFO: SRF file: After scaling time series: dt*sum(par)=%e [m]\n", slip_sum);
        }
      }
      //done scaling        
      
      npar = nt1dim+2;
      nipar = 1;
  
      // printf("Read discrete time series: tinit=%g, dt=%g, nt1=%i\n", tinit, dt, nt1);
      // for (int i=0; i<nt1+1; i++)
      //   printf("Sv1[%i]=%g\n", i+1, par[i+1]);
  
  // convert lat, lon, depth to (x,y,z)
      ew->computeCartesianCoord(x, y, lon, lat);
  // convert depth in [km] to [m]
      z = dep * 1e3;
  
  // convert strike, dip, rake to Mij
      float_sw4 radconv = M_PI / 180.;
      float_sw4 S, D, R;
      stk -= mGeoAz; // subtract off the grid azimuth
      S = stk*radconv; D = dip*radconv; R = rake*radconv;
    
      mxx = -1.0 * ( sin(D) * cos(R) * sin (2*S) + sin(2*D) * sin(R) * sin(S)*sin(S) );
      myy =        ( sin(D) * cos(R) * sin (2*S) - sin(2*D) * sin(R) * cos(S)*cos(S) );
      mzz = -1.0 * ( mxx + myy );	
      mxy =        ( sin(D) * cos(R) * cos (2*S) + 0.5 * sin(2*D) * sin(R) * sin(2*S) );
      mxz = -1.0 * ( cos(D) * cos(R) * cos (S)   + cos(2*D) * sin(R) * sin(S) );
      myz = -1.0 * ( cos(D) * cos(R) * sin (S)   - cos(2*D) * sin(R) * cos(S) );
  
  // scale (note that the shear modulus is not yet available. Also note that we convert [cm] to [m])
      m0 = area*1e-4 * slip1*1e-2;
    
      mxx *= m0;
      mxy *= m0;
      mxz *= m0;
      myy *= m0;
      myz *= m0;
      mzz *= m0;
  
  // before creating the source, make sure (x,y,z) is inside the computational domain
  
  // only check the z>zmin when we have topography. For a flat free surface, we will remove sources too 
  // close or above the surface in the call to mGlobalUniqueSources[i]->correct_Z_level()
  
      if (x < xmin || x > m_global_xmax || y < ymin || y > m_global_ymax || z < zmin || z > m_global_zmax)
      {
        stringstream sourceposerr;
        sourceposerr << endl
                     << "***************************************************" << endl
                     << " ERROR:  Source positioned outside grid!  " << endl
                     << endl
                     << " Source from rupture file @" << endl
                     << "  x=" << x << " y=" << y << " z=" << z << endl 
                     << "  lat=" << lat << " lon=" << lon << " dep=" << dep << endl 
                     << endl;
          
        if ( x < xmin )
          sourceposerr << " x is " << xmin - x << 
            " meters away from min x (" << xmin << ")" << endl;
        else if ( x > m_global_xmax)
          sourceposerr << " x is " << x - m_global_xmax << 
            " meters away from max x (" << m_global_xmax << ")" << endl;
        if ( y < ymin )
          sourceposerr << " y is " << ymin - y << 
            " meters away from min y (" << ymin << ")" << endl;
        else if ( y > m_global_ymax)
          sourceposerr << " y is " << y - m_global_ymax << 
            " meters away from max y (" << m_global_ymax << ")" << endl;
        if ( z < zmin )
          sourceposerr << " z is " << zmin - z << 
            " meters away from min z (" << zmin << ")" << endl;
        else if ( z > m_global_zmax)
          sourceposerr << " z is " << z - m_global_zmax << 
            " meters away from max z (" << m_global_zmax << ")" << endl;
        sourceposerr << "***************************************************" << endl;
        if (world_rank == 0)
          cout << sourceposerr.str();
      }
      else if( !skip_zero_slip_point )
      {
        sourcePtr = new Source(ew, freq, t0, x, y, z, mxx, mxy, mxz, myy, myz, mzz,
                               tDep, formstring, topodepth, ncyc, par, npar, ipar, nipar, true ); // true is correctStrengthForMu
  
        if (sourcePtr->ignore())
        {
          delete sourcePtr;
        }
        else
        {
          a_GlobalUniqueSources[event].push_back(sourcePtr);
          nSources++;
        }
      }
  
      // deallocate temporary arrays...
      delete[] par;
      delete[] ipar;
  
    } // end if nt1 >0
  
    // read past discrete time series for u2
    if (nt2>0)
    {
      nu2++;
      double dum;
      if (world_rank == 0)
        printf("WARNING nt2=%i > 0 will be ignored\n", nt2);
    } // end if nt2 > 0
  
    // read past discrete time series for u3
    if (nt3>0)
    {
      nu3++;
      double dum;
      if (world_rank == 0)
        printf("WARNING nt3=%i > 0 will be ignored\n", nt3);
    } // end if nt3 > 0
    
  } // end for all sources
  if (world_rank == 0)
    printf("Read npts=%i, made %i point moment tensor sources, nu1=%i, nu2=%i, nu3=%i\n", npts, nSources, nu1, nu2, nu3);
  if (world_rank == 0 && nskip_zero_slip > 0)
    printf("Skipped %i rupture points with zero slip-velocity integral in u1.\n", nskip_zero_slip);

  etime = MPI_Wtime();
  if (is_debug && world_rank == 0) 
      printf("Create source takes %.2f seconds\n", etime-stime);


  H5Tclose(ctype);
  H5Tclose(dtype);
  free(point_data);
  free(sr_data);

  return;
}

#endif // USE_HDF5
#endif // READHDF5_C
