from numpy import *
from random import uniform
import tqdm
from random import  uniform
#import sys
import time
from multiprocessing import Pool
from numpy import linspace,sqrt,zeros
from tqdm import tqdm_notebook
from random import  uniform
from multiprocessing import Pool
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import AutoMinorLocator
import matplotlib.pyplot as plt
plt.matplotlib.rcParams.update({'font.size': 30,'legend.fontsize':25})#,'font.family': 'serif'})
plt.rc('axes', linewidth=3)
#plt.rc('font.size',25)


def make_outputfile(fig,tt='test',res=200):
    fname_eps = 'Fig_{t1}.eps'.format(t1=tt)
    fname_png = 'Fig_{t1}.png'.format(t1=tt)
    fname_pdf = 'Fig_{t1}.pdf'.format(t1=tt)
    fname_svg = 'Fig_{t1}.svg'.format(t1=tt)
    fig.savefig(fname_eps,bbox_inches='tight', transparent=True)
    fig.savefig(fname_png,dpi=res,bbox_inches='tight', transparent=True)
    fig.savefig(fname_pdf,dpi=res,bbox_inches='tight', transparent=True)
#    fig.savefig(fname_svg,dpi=res,bbox_inches='tight', transparent=True)
    print("Files '{t1}', '{t3}', '{t2} and {t5}' created, where the latter two have dpi={t4}".format(t1=fname_eps,t2=fname_png,t3=fname_pdf,t5=fname_svg,t4=str(res)))
    
    
def read_input(inputfile):
    """
    Reading input data from bouler.in
    """
    from numpy import append
#    inputfile = 'boulder.in'
    with open(inputfile) as f:
        lines = f.readlines()
        dummy =[]
        for line in lines:
            if line.split()[0] != 'slope:' and line.split()[0] != 'delta:':
                dummy.append(float(line.split()[-1]))
    return dummy


def read_input_sldl(inputfile):
    with open(inputfile) as f:
        lines = f.readlines()
        slope_ll = []
        delta_ll = []
        for line in lines:
            if line.split()[0] == 'slope:':
                for i in range(1,len(line.split(',')[:])-1):
                    slope_ll.append(float((line.split(',')[i])))
            if line.split()[0] == 'delta:':
                for i in range(1,len(line.split(',')[:])-1):
                    delta_ll.append(float(line.split(',')[i]))
    delta_ll = array(delta_ll)
    slope_ll = array(slope_ll)
    return slope_ll, delta_ll

def read_inputparastudy(parafile,parameter_l):
    print('in function:',parameter_l)
    check_text = '{t1}:'.format(t1=parameter_l)
    output_array = []
    if parameter_l != 'Range':
        print('wrong loop')
        with open(parafile) as f:
            lines = f.readlines()
            for line in lines:
                print(line)
                if line.split()[0] == check_text:
                    for i in range(1,len(line.split(',')[:])-1):
                        output_array.append(float(line.split(',')[i]))
        print(output_array)   
        return output_array 
     
    if parameter_l == 'Range':
        print('in loop')
        output_array_1 = []
        output_array_2 = []
        output_array_3 = []
        with open(parafile) as f:
            lines = f.readlines()
            for line in lines:
                print(line)
                if line.split()[0] == 'Density:':
                    for i in range(1,len(line.split(',')[:])-1):
                        output_array_1.append(float(line.split(',')[i]))
                if line.split()[0] == 'Size:':
                    for i in range(1,len(line.split(',')[:])-1):
                        output_array_2.append(float(line.split(',')[i]))
                if line.split()[0] == 'MC:':
                    for i in range(1,len(line.split(',')[:])-1):
                        output_array_3.append(float(line.split(',')[i]))
        print(output_array_1,output_array_2,output_array_3)
        return output_array_3[0],output_array_1[0],output_array_1[1],output_array_2[0],output_array_2[1]
def read_outputparastudy(parafile,parameter_l):
    check_text = '{t1}:'.format(t1=parameter_l)
    output_array = []
    with open(parafile) as f:
        lines = f.readlines()
        for line in lines:
            if line.split()[0] == check_text:
                for i in range(1,len(line.split(',')[:])-1):
                    output_array.append(str(line.split(',')[i]))
   
    print(output_array)   
    return output_array 




def read_run_para(parafile):
    with open(parafile) as f:
        lines = f.readlines()
        for line in lines:
            print(line)
            if line.split()[0] == 'Simulation-Parameter-File:':
                para_filePath = str(line.split()[-1])
            if line.split()[0] == 'Platform:':
                platformName = str(line.split()[-1])
            if line.split()[0] == 'Target-Folder:':
                target_folder = str(line.split()[-1])
    return platformName, para_filePath,target_folder

def data_folder_checking(path_l):
    import os
    import glob
    testing = os.path.exists(path_l)
    print(testing)
    if testing == False:
        os.mkdir(path_l)
    if testing == True:
        num_files = len(glob.glob(path_l+'*'))
        if num_files > 0:
            files = glob.glob(path_l+'*')
            for f in files:
                    os.remove(f)



def round_n(num, sig_figs):
    """Round to specified number of sigfigs.     
    """
    if num != 0:
        return round(num, -int(math.floor(math.log10(abs(num))) - (sig_figs - 1)))
    else:
        return 0  # Can't take the log of 0

def motion_space_p(fr_l,fd_l,dis_l,ca):
    """
    Analyzes the boulder simulation results and checks if for any given
    realization, a boulder dislodged or not (for both, sliding and rotation).
    """
    len_mc_l = len(fr_l)
    fr_dis_l = []
    fd_dis_l = []
    fr_ndis_l = []
    fd_ndis_l = []
    for i_mc in range(len_mc_l):
        if max(dis_l[:,i_mc]) >= ca:
            fr_dis_l.append(fr_l[i_mc])
            fd_dis_l.append(fd_l[i_mc])
        if max(dis_l[:,i_mc]) < ca:
            fr_ndis_l.append(fr_l[i_mc])
            fd_ndis_l.append(fd_l[i_mc])
    fr_dis_l = array(fr_dis_l)
    fd_dis_l = array(fd_dis_l)
    fr_ndis_l = array(fr_ndis_l)
    fd_ndis_l = array(fd_ndis_l)
    return fr_dis_l,fd_dis_l,fr_ndis_l,fd_ndis_l



def func(index_total):
    i = int(floor(index_total/len(dl_changeg)))
    j = index_total % len(dl_changeg)
    filename = 'out_'+det_runname(round_n(sl_changeg[i],3),round_n(dl_changeg[j],3))
    #fl1 = './data/'+filename+'_conditions.dat'
    #fl2 = './data/'+filename+'_rotation.dat'
    fl1 = path_g+filename+'_conditions.dat'
    fl2 = path_g+filename+'_rotation.dat'
    #fl1= 'out_sln1000_dlp0150_conditions.dat'
    #fl2= 'out_sln1000_dlp0150_rotation.dat'

    #print(fl1,path.exists(fl1))
   # print(fl2,path.exists(fl2))
   # print(m,sl_i,j,filename)
    data_dis = loadtxt(fl2)#,delimiter=",")
    condi = loadtxt(fl1)#,delimiter=",")
    fd = condi[:,1]
    fr = condi[:,0]
    time = data_dis[:,0]
    critical_angle = 0.5 * pi
    critical_dista = 1.0
    fr_dis ,fd_dis,fr_ndis ,fd_ndis = \
        motion_space_p(fr,fd,data_dis[:,1:],critical_angle)
    return index_total,fr_dis,fd_dis,fr_ndis,fd_ndis

def read_model_data(path_l='./',environ='conv'):
    """
    Reading data saved in the given path. The program assumes that the input
    file for the boulder simulations are inside the directory
    """
    import platform
    import tqdm
    import glob
    import sys
    platform.platform()
    platform.system()
    input_file_name_l = glob.glob(path_l+'bould*.in')
    print(input_file_name_l)
    if len(input_file_name_l) == 0:
        print('No input file present!')
        return 0.0,0.0,0.0,0.0,0.0,0.0
    if len(input_file_name_l) > 1:
        print('Too many input files:')
        return
    if len(input_file_name_l) == 1:
        d1,d1,d1,MC_set,d1,d1,d1,d1,d1,d1,d1,d1,d1,d1,d1,d1,d1,\
                d1,d1,d1,d1=read_input(input_file_name_l[0])
        mc = int(MC_set)
#d1,d1,d1,MC_set,d2=read_input('./boulder.in')
        sl_changel, dl_changel = read_input_sldl(input_file_name_l[0])
        nx = len(dl_changel)*len(sl_changel)
        c_fu_disl =  zeros([nx,mc])
        c_fd_disl =  zeros([nx,mc])
        c_fu_ndisl =  zeros([nx,mc])
        c_fd_ndisl =  zeros([nx,mc])
        if platform.system() == 'Linux':
            print('\t\t \033[1mLinux\033[0m')
            from multiprocessing import Pool,cpu_count
            from functools import partial
            inputs = range(nx)
            global sl_ig,path_g
            global sl_changeg,dl_changeg
            sl_changeg = sl_changel
            dl_changeg = dl_changel
            n_proc = cpu_count()
            list_index = []
            mm=-1
            path_g = path_l
            print('Reading data in parallel using {t1} CPUs (max amount on machine).'.format(t1=n_proc))
            pool = Pool(processes = n_proc)
            if environ == 'noteb':
                for res in tqdm.tqdm_notebook(pool.imap(func, inputs), total=nx,desc = 'Processing'):
                    c_fu_disl[res[0],0:len(res[1])] = res[1][:]
                    c_fd_disl[res[0],0:len(res[2])] = res[2][:]
                    c_fu_ndisl[res[0],0:len(res[3])] = res[3][:]
                    c_fd_ndisl[res[0],0:len(res[4])] = res[4][:]


            else:
                for res in tqdm.tqdm(pool.imap(func, inputs), total=nx,desc='Processing'):
                    c_fu_disl[res[0],0:len(res[1])] = res[1][:]
                    c_fd_disl[res[0],0:len(res[2])] = res[2][:]
                    c_fu_ndisl[res[0],0:len(res[3])] = res[3][:]
                    c_fd_ndisl[res[0],0:len(res[4])] = res[4][:]
            pool.close()
            pool.terminate()
            return sl_changel,dl_changel,c_fu_disl,c_fd_disl,c_fu_ndisl,c_fd_ndisl
        if platform.system() == 'Darwin':
            print(platform.sys())
            m=-1
            for i in range(nx):
                sl_i = int(floor(i/16))
                j= i%16
                filename = 'out_'+det_runname(round_n(sl_changel[sl_i],3),round_n(dl_changel[j],3))
                #fl1 = './data/'+filename+'_conditions.dat'
                #fl2 = './data/'+filename+'_rotation.dat'
                fl1 = filename+'_conditions.dat'
                fl2 = filename+'_rotation.dat'
                data_dis = loadtxt(fl2)#,delimiter=",")
                condi = loadtxt(fl1)#,delimiter=",")
                fd = condi[:,1]
                fr = condi[:,0]
                time = data_dis[:,0]
                critical_angle = 0.5 * pi
                critical_dista = 1.0
                fr_dis ,fd_dis,fr_ndis ,fd_ndis = \
                    motion_space_p(fr,fd,data_dis[:,1:],critical_angle)
                c_fu_disl[i,0:len(fr_dis)] = fr_dis[:]
                c_fd_disl[i,0:len(fd_dis)] = fd_dis[:]
                c_fu_ndisl[i,0:len(fr_ndis)] = fr_ndis[:]
                c_fd_ndisl[i,0:len(fd_ndis)] = fd_ndis[:]
            return sl_changel,dl_changel,c_fu_disl,c_fd_disl,c_fu_ndisl,c_fd_ndisl


def read_input_RandomBlocks():
    """
    Reading input data from bouler.in
    """
    from numpy import append
    inputfile = 'boulder_random.in'
    with open(inputfile) as f:
        lines = f.readlines()
        dummy =[]
        for line in lines:
            dummy.append(float(line.split()[-1]))
    return dummy



def read_data(fname):
    from numpy import loadtxt
    array1 = loadtxt(fname)
    return array1


def replace_line(file_name, line_num, text):
    """
    Replace line in respective file. Please note that the line
    numbering starts at zero.
    """
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()
    
def round_n(num, sig_figs):
    """Round to specified number of sigfigs.
    """
    if num != 0:
        return round(num, -int(math.floor(math.log10(abs(num))) - (sig_figs - 1)))
    else:
        return 0  # Can't take the log of 0    

def det_runname(slope_l,delta_l):
    if slope_l < 0.0:
        test_sl = 'n{t1:04d}'.format(t1=abs(int(slope_l*100)))
    else:
        test_sl = 'p{t1:04d}'.format(t1=abs(int(slope_l*100)))
    if delta_l < 0.0:
        test_dl = 'n{t1:04d}'.format(t1=abs(int(delta_l*100)))
    else:
        test_dl = 'p{t1:04d}'.format(t1=abs(int(delta_l*100)))
    eqname = 'sl{t1}_dl{t2}'.format(t1=test_sl,t2=test_dl)
    return eqname    
    
def debugplots(x2,x3,x4,z2,z3,z4,xv,zv,xv1,zv1,xv2,zv2,xv3,zv3,xv4,zv4):
    """
    Check grid cells retained in blockmoments is correct, for debugging
    """
    import matplotlib.pyplot as plt
    print('plotting results')
    fig, (ax1,ax2) = plt.subplots(nrows=1,ncols=2,figsize=(22,10))
    ax1.plot([0.,x2,x3,x4,0.],[0.,z2,z3,z4,0.],'k')
    ax1.scatter(xv,zv,marker='o',color='gray',s=10,alpha=0.2)
    ax1.scatter(xv1,zv1,marker='o',color='blue',s=7,alpha=0.2)
    ax1.scatter(xv2,zv2,marker='o',color='red',s=5,alpha=0.2)
    ax1.set_xlabel('x [m]');
    ax1.set_ylabel('z [m]');
    ax1.axis('equal')
    ax2.plot([0.,x2,x3,x4,0.],[0.,z2,z3,z4,0.],'k')
    ax2.scatter(xv3,zv3,marker='o',color='green',s=10,alpha=0.2)
    ax2.scatter(xv4,zv4,marker='o',color='black',s=7,alpha=0.2)
    ax2.set_xlabel('x [m]');
    ax2.set_ylabel('z [m]');
    ax2.axis('equal')
    plt.show()

   
def save_md_data(filename,arr):
    savetxt(filename,arr,fmt='%5.4e',delimiter=' ')
#    with open(filename,'w') as f:
#        for row in arr:
#            row.tofile(f,format='%5.4e',sep=',')
#            f.write('\n')    
    
def save_data(arr,file_name):
    from numpy import savetxt
    outfile=open(file_name,'w')
    savetxt(outfile, arr)
    outfile.close()

   
def plot_data(fname,aaa,bbb,ccc,n_froude,n_eta,rhos_min,rhos_max):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from numpy import shape,vander,linalg,poly1d
    ar = read_data(fname)
    degree = 8
    A = vander(ar[:,0], degree)
    (coeffs, residuals, rank, sing_vals) = linalg.lstsq(A, ar[:,1])
    f = poly1d(coeffs)
    y_est = f(ar[:,0])
    plt.scatter(ar[:,0],ar[:,1],marker='x',s=20,c='k',lw=0.5,label='Data')
    #plt.plot(ar[:,0],smooth(ar[:,1],50),'r-')
    plt.plot(ar[:,0],y_est,'r-',lw=1,label = 'Fit')
    plt.legend()
    plt.xlabel('Froude Number')
    plt.ylabel('Flow Depth [m]')
    plt.xlim(0.5,2.5)
    plt.ylim(0,20)
    plt.title("Flow depth to flip over a boulder ({t1} x {t2} x {t3}) as a f(Fr)".format(t1 =aaa, t2 = bbb, t3 = ccc))
    plt.text(0.55,19.2,"This diagram contains {t1} individual runs".format(t1=n_froude*n_eta))
    #    plt.savefig('test_bahamas.png',dpi=501, bbox_inches="tight",transparent=True)
    plt.savefig('test_bahamas.png',dpi=501, bbox_inches="tight")         
