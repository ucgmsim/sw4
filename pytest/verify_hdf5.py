import numpy as np
import h5py
import struct
import os

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

def read_sac_usgs(fname):
    x = []
    y = []
    z = []
    t = []
    sta = open(fname,'r')
    lines = sta.readlines()
    skip = 0
    flag = 0
    count = 0

    for line in lines:
        if line[0] == '#':
            skip += 1
            continue
        elif flag == 0:
            n = len(lines) - skip
            flag = 1
            x = np.zeros(n)
            y = np.zeros(n)
            z = np.zeros(n)
            t = np.zeros(n)
            
        t[count] = np.float32(line.split()[0])
        x[count] = np.float32(line.split()[1])
        y[count] = np.float32(line.split()[2])
        z[count] = np.float32(line.split()[3])
        count += 1
    sta.close()
    return (np.array(t), np.array(x), np.array(y), np.array(z))

def read_sac_hdf5(fname, staname):
    fid = h5py.File(fname, 'r')
    grp = fid[staname]
    npts = grp['NPTS'][0]
    nsew = grp['ISNSEW'][0]
    delta = fid['DELTA'][0]
    
    if nsew == 1:
        data_names = ['EW', 'NS', 'UP']
        has_nsew = 1
    else:
        data_names = ['X', 'Y', 'Z']
        
    x=np.array(grp[data_names[0]][0:npts])
    y=np.array(grp[data_names[1]][0:npts])
    z=np.array(grp[data_names[2]][0:npts])
    
    if npts <= 0 or delta <=0:
        print('Error with station data, npts = ', npts, ', delta = ', delta)
    origintime = fid['ORIGINTIME'][0]
    t = np.linspace(0, delta*(npts-1), npts)
        
    fid.close()
    return (t, x, y, z)

def read_sw4img(fname):
    img = open(fname,'rb')
    prec = struct.unpack('i', img.read(4))[0]
    npatch = struct.unpack('i', img.read(4))[0]
    t = struct.unpack('d', img.read(8))[0]
    plane = struct.unpack('i', img.read(4))[0]
    xi = struct.unpack('d', img.read(8))[0]
    mode = struct.unpack('i', img.read(4))[0]
    grid_info = struct.unpack('i', img.read(4))[0]
    # Skip create time
    img.read(25)
    
    #print("prec=%d, npatch=%d, t=%f, plane=%d, xi=%f, mode=%d, grid_info=%d" % \
    #      (prec, npatch, t, plane, xi, mode, grid_info))
    
    # Grid size info
    h = np.zeros(npatch, dtype=np.float64)
    zmin = np.zeros(npatch, dtype=np.float64)
    i = np.zeros(npatch, dtype=int)
    ni = np.zeros(npatch, dtype=int)
    j = np.zeros(npatch, dtype=int)
    nj = np.zeros(npatch, dtype=int)
    nelem = 0
    for u in range (0, npatch):
        h[u] = struct.unpack('d', img.read(8))[0]
        zmin[u] = struct.unpack('d', img.read(8))[0]
        i[u] = struct.unpack('i', img.read(4))[0]
        ni[u] = struct.unpack('i', img.read(4))[0]
        j[u] = struct.unpack('i', img.read(4))[0]
        nj[u] = struct.unpack('i', img.read(4))[0]
        nelem += ni[u] * nj[u]
        
        #print("patch %d: h=%f, zmin=%f, i=%d, ni=%d, j=%d, nj=%d" % \
        #      (u, h[u], zmin[u], i[u], ni[u], j[u], nj[u]))

    # Read all patch data
    if prec == 4:
        pdata = struct.unpack(str(nelem)+'f', img.read(prec*nelem))
    elif prec == 8:
        pdata = struct.unpack(str(nelem)+'d', img.read(prec*nelem))
    #print("patch %d: min=%e, max=%e" % (u, np.min(pdata), np.max(pdata)))
    img.close()
    return np.array(pdata)

def read_sw4img_hdf5(fname):
        
    image_h5      = h5py.File(fname, 'r')
    npatch        = image_h5['npatch'][0]
    time          = image_h5['time'][0]
    plane         = image_h5['plane'][0]
    coordinate    = image_h5['coordinate'][0]
    grid_info     = image_h5['gridinfo'][0]
    mode          = image_h5['mode'][0]
    creation_time = str(image_h5.attrs['creationtime'], 'utf-8')

    grid_sizes    = image_h5['grid_size']
    zmins         = image_h5['zmin']
    nis           = image_h5['ni']
    njs           = image_h5['nj']
    
    readz = False
    has_grid = False
    if 'grid' in image_h5.keys():
        has_grid = True
    
    pdata = image_h5['patches']

    return np.array(pdata)


def read_essi(fname):
    essi = h5py.File(fname, 'r')
    data0 = essi["vel_0 ijk layout"][:]
    data1 = essi["vel_1 ijk layout"][:]
    data2 = essi["vel_2 ijk layout"][:]
    
    essi.close()
    return data0, data1, data2


def verify_supergrid_metadata(fname, sta_names, width_m, width_gp,
                              require_present=True):
    """Check the supergrid metadata SW4 writes into a station file.

    SGDEPTH is 0 for a station in the interior, which is a positive statement
    rather than a default: the datasets are created with H5D_FILL_TIME_NEVER,
    so an unwritten one reads back as garbage and this catches that mode.
    width_m/width_gp are the expected file-level SGWIDTH/SGWIDTHGP.
    require_present=False tolerates their absence, for a restart run whose
    output file was created by an earlier invocation.
    """
    fid = h5py.File(fname, 'r')
    ok = True
    try:
        have = 'SGWIDTH' in fid and 'SGWIDTHGP' in fid
        if not have:
            if require_present:
                print('Supergrid metadata: file-level SGWIDTH/SGWIDTHGP missing')
                return False
            return True
        got_m = float(np.array(fid['SGWIDTH']).squeeze())
        got_gp = float(np.array(fid['SGWIDTHGP']).squeeze())
        if abs(got_m - width_m) > 1e-6 * max(1.0, abs(width_m)):
            print('Supergrid metadata: SGWIDTH = %g, expected %g' % (got_m, width_m))
            ok = False
        if abs(got_gp - width_gp) > 1e-6 * max(1.0, abs(width_gp)):
            print('Supergrid metadata: SGWIDTHGP = %g, expected %g' % (got_gp, width_gp))
            ok = False
        for sta in sta_names:
            grp = fid[sta]
            for name in ('SGDEPTH', 'SGDEPTHGP'):
                if name not in grp:
                    print('Supergrid metadata: [%s/%s] missing' % (sta, name))
                    ok = False
                    continue
                val = float(np.array(grp[name]).squeeze())
                if val != 0.0:
                    print('Supergrid metadata: [%s/%s] = %g, expected 0 (all ten '
                          'loh1 stations are in the interior)' % (sta, name, val))
                    ok = False
    finally:
        fid.close()
    return ok

def verify(pytest_dir, tolerance):
    ref_dir = pytest_dir + '/hdf5/loh1-h100-mr-1/'
    hdf5_dir = os.getcwd() + '/loh1-h100-mr-1-hdf5/'
    verify = True
    nsta = 0
    for i in range(1,11):
        sta_name = 'sta%02d' % i
        #print(sta_name)
        usgs_fname = ref_dir + sta_name + '.txt'
        hdf5_fname = hdf5_dir + 'sta.h5'
        usgs_t, usgs_x, usgs_y, usgs_z = read_sac_usgs(usgs_fname)
        hdf5_t, hdf5_x, hdf5_y, hdf5_z = read_sac_hdf5(hdf5_fname, sta_name)
        if np.max(hdf5_t-usgs_t) > tolerance or np.min(hdf5_t-usgs_t) < -tolerance:
            verify = False
            print ("Station [%s] time data not match!" % sta_name)
            return False
        if np.max(hdf5_x-usgs_x) > tolerance or np.min(hdf5_x-usgs_x) < -tolerance:
            verify = False
            print ("Station [%s] x data not match!" % sta_name)
            print(hdf5_x-usgs_x)
            return False
        if np.max(hdf5_y-usgs_y) > tolerance or np.min(hdf5_y-usgs_y) < -tolerance:
            verify = False
            print ("Station [%s] y data not match!" % sta_name)
            return False
        if np.max(hdf5_z-usgs_z) > tolerance or np.min(hdf5_z-usgs_z) < -tolerance:
            verify = False
            print ("Station [%s] z data not match!" % sta_name)  
            return False
        nsta += 1

    # if verify == 1:
    #     print ('All %d stations data match!' % nsta)

    # h=200, supergrid gp=30 -> a 6000 m layer, 30 grid points.
    if not verify_supergrid_metadata(hdf5_fname,
                                     ['sta%02d' % i for i in range(1, 11)],
                                     6000.0, 30.0):
        return False

    nimg = 0
    for filename in os.listdir(ref_dir):
        if filename.endswith(".sw4img"):
            nimg += 1
            #print(filename)
            sw4img_file = ref_dir + filename
            h5img_file  = hdf5_dir + filename + '.h5'
            sw4_pdata = read_sw4img(sw4img_file)
            h5_pdata = read_sw4img_hdf5(h5img_file)
            if len(h5_pdata) != len(h5_pdata):
                print("Image sizes are diferent! %d/%d" % (len(sw4_pdata), len(h5_pdata)))
                verify = False
                return False
            else:
                if np.max(h5_pdata-sw4_pdata) > tolerance or np.min(h5_pdata-sw4_pdata) < -tolerance:
                    print ("Image data [%s] does not match!" % sw4img_file)        
                    verify = False
                    return False
                    
    # if verify == 1:
    #     print ('All %d images data match!' % nimg)

    essi_fname = hdf5_dir + 'essioutput.ssi'
    ref_essi_fname = ref_dir + 'essioutput.cycle=000.essi'
    data0, data1, data2 = read_essi(essi_fname)
    ref_data0, ref_data1, ref_data2 = read_essi(ref_essi_fname)
    sum0 = np.sum(data0)
    sum1 = np.sum(data1)
    sum2 = np.sum(data2)
    min0 = np.min(data0)
    min1 = np.min(data1)
    min2 = np.min(data2)
    max0 = np.max(data0)
    max1 = np.max(data1)
    max2 = np.max(data2)

    ref_sum0 = np.sum(ref_data0)
    ref_sum1 = np.sum(ref_data1)
    ref_sum2 = np.sum(ref_data2)
    ref_min0 = np.min(ref_data0)
    ref_min1 = np.min(ref_data1)
    ref_min2 = np.min(ref_data2)
    ref_max0 = np.max(ref_data0)
    ref_max1 = np.max(ref_data1)
    ref_max2 = np.max(ref_data2)
    # larger tolerance for sum data
    if np.absolute(sum0-ref_sum0) > tolerance*1e4 or np.absolute(sum1-ref_sum1) > tolerance*1e4 or np.absolute(sum2-ref_sum2) > tolerance*1e4:
        print ("ESSI data sum not match!", sum0, sum1, sum2)
        return False

    if np.absolute(min0-ref_min0) > tolerance or np.absolute(min1-ref_min1) > tolerance or np.absolute(min2-ref_min2) > tolerance:
        print ("ESSI data min not match!", min0, min1, min2)
        return False

    if np.absolute(max0-ref_max0) > tolerance or np.absolute(max1-ref_max1) > tolerance or np.absolute(max2-ref_max2) > tolerance:
        print ("ESSI data max not match!", max0, max1, max2)
        return False

    return verify

def verify_sac_image(pytest_dir, tolerance):
    ref_dir = pytest_dir + '/hdf5/loh1-h100-mr-1/'
    hdf5_dir = os.getcwd() + '/loh1-h100-mr-1-hdf5/'
    verify = True
    nsta = 0
    for i in range(1,11):
        sta_name = 'sta%02d' % i
        #print(sta_name)
        usgs_fname = ref_dir + sta_name + '.txt'
        hdf5_fname = hdf5_dir + 'sta.h5'
        usgs_t, usgs_x, usgs_y, usgs_z = read_sac_usgs(usgs_fname)
        hdf5_t, hdf5_x, hdf5_y, hdf5_z = read_sac_hdf5(hdf5_fname, sta_name)
        if np.max(hdf5_t-usgs_t) > tolerance or np.min(hdf5_t-usgs_t) < -tolerance:
            verify = False
            print ("Station [%s] time data not match!" % sta_name)
        if np.max(hdf5_x-usgs_x) > tolerance or np.min(hdf5_x-usgs_x) < -tolerance:
            verify = False
            print ("Station [%s] x data not match!" % sta_name)
            print(hdf5_x-usgs_x)
        if np.max(hdf5_y-usgs_y) > tolerance or np.min(hdf5_y-usgs_y) < -tolerance:
            verify = False
            print ("Station [%s] y data not match!" % sta_name)
        if np.max(hdf5_z-usgs_z) > tolerance or np.min(hdf5_z-usgs_z) < -tolerance:
            verify = False
            print ("Station [%s] z data not match!" % sta_name)  
        nsta += 1

    # if verify == 1:
    #     print ('All %d stations data match!' % nsta)

    # Same geometry; tolerate absent metadata because the restart case reuses a
    # station file created by an earlier invocation.
    if not verify_supergrid_metadata(hdf5_fname,
                                     ['sta%02d' % i for i in range(1, 11)],
                                     6000.0, 30.0, require_present=False):
        return False

    nimg = 0
    for filename in os.listdir(ref_dir):
        if filename.endswith(".sw4img"):
            nimg += 1
            #print(filename)
            sw4img_file = ref_dir + filename
            h5img_file  = hdf5_dir + filename + '.h5'
            sw4_pdata = read_sw4img(sw4img_file)
            h5_pdata = read_sw4img_hdf5(h5img_file)
            if len(h5_pdata) != len(h5_pdata):
                print("Image sizes are diferent! %d/%d" % (len(sw4_pdata), len(h5_pdata)))
                verify = False
            else:
                if np.max(h5_pdata-sw4_pdata) > tolerance or np.min(h5_pdata-sw4_pdata) < -tolerance:
                    print ("Image data [%s] does not match!" % sw4img_file)        
                    verify = False
            if verify == False:
                break
                    
    # if verify == 1:
    #     print ('All %d images data match!' % nimg)

    return verify
