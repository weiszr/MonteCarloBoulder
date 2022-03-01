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
from helpers import make_outputfile
plt.matplotlib.rcParams.update({'font.size': 30,'legend.fontsize':25})#,'font.family': 'serif'})
plt.rc('axes', linewidth=3)
#plt.rc('font.size',25)
from helpers import *

def motion_space_c(fr_l,fd_l,sli_l,dis_l,cd,ca):
    len_mc_l = len(fr_l)
    fr_slidis_l = []
    fd_slidis_l = []
    fr_slindis_l = []
    fd_slindis_l = []
    fr_nslindis_l = []
    fd_nslindis_l = []
    for i_mc in range(len_mc_l):
        if max(sli_l[:,i_mc]) >= cd and max(dis_l[:,i_mc]) >= ca:
        #if max(dis_l[:,i_mc]) >= ca:
            fr_slidis_l.append(fr_l[i_mc])
            fd_slidis_l.append(fd_l[i_mc])
        if max(sli_l[:,i_mc]) < cd and max(dis_l[:,i_mc]) < ca:
            fr_nslindis_l.append(fr_l[i_mc])
            fd_nslindis_l.append(fd_l[i_mc])
        if max(sli_l[:,i_mc]) >= cd and max(dis_l[:,i_mc]) < ca:
            fr_slindis_l.append(fr_l[i_mc])
            fd_slindis_l.append(fd_l[i_mc])
    fr_slidis_l = array(fr_slidis_l)
    fd_slidis_l = array(fd_slidis_l)
    fr_nslindis_l = array(fr_nslindis_l)
    fd_nslindis_l = array(fd_nslindis_l)
    fr_slindis_l = array(fr_slindis_l)
    fd_slindis_l = array(fd_slindis_l)
    print(shape(fr_slidis_l))
    print(shape(fr_nslindis_l))
    print(shape(fr_slindis_l))
    return fr_slidis_l,fd_slidis_l,fr_nslindis_l,fd_nslindis_l,fr_slindis_l,fd_slindis_l

def motion_space_p(fr_l,fd_l,dis_l,ca):
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



def find_boundary_impr(arr_x_l,arr_y_l, arr_z_l1,arr_z_l2,res):    
    nx = len(arr_x_l)
    fac_res =10
    radius = 1.0
    new_arr_x_l = linspace(min(arr_x_l),max(arr_x_l),res)
    new_arr_y_l = linspace(min(arr_y_l),max(arr_y_l),fac_res*res)
    val1 = zeros([res,fac_res*res])
    val2 = zeros([res,fac_res*res])
    for i_res in range(res):
        for j_res in range(fac_res*res):
            dist = []
            index_list = []            
            for i_arr in range(nx):
#                print(i_res,j_res,i_arr)
                sq_dist = (new_arr_x_l[i_res] - arr_x_l[i_arr])**2.0 + \
                                (new_arr_y_l[j_res] - arr_y_l[i_arr])**2.0
                    # print(sq_dist)
                if sq_dist <= radius**2.0:
                    index_list.append(i_arr)
                    dist.append(sq_dist)
            dist = array(dist)
#            print(dist)
            index_temp = argmin(dist)
            print(index_temp,nx,shape(arr_z_l1))
            val1[i_res,j_res] = max(arr_z_l1[:,index_temp])
    for i_res in range(res):
        for j_res in range(fac_res*res):
            dist = []
            index_list = []
            for i_arr in range(nx):
                sq_dist = (new_arr_x_l[i_res] - arr_x_l[i_arr])**2.0 + \
                                (new_arr_y_l[j_res] - arr_y_l[i_arr])**2.0
                    # print(sq_dist)
                if sq_dist <= radius**2.0:
                    index_list.append(i_arr)
                    dist.append(sq_dist)
            dist = array(dist)
            index_temp = argmin(dist)
            val2[i_res,j_res] = max(arr_z_l2[:,index_temp])
#    for i_res in range(res):
#        for j_res in range(fac_res*res):
#            if val1[i_res,j_res] < 0.5 *pi and   
    print(val1)      
    return new_arr_x_l,new_arr_y_l,val1,val2
    
            # dens_output = zeros(res)
    # if type == 'max':
    #     for i_res in range(res):
    #         m=0
    #         for j_res in reversed(range(fac_res*res)):
    # 
    #                 len_mc_l = len(fr_l)
    # fr_slidis_l = []
    # fd_slidis_l = []
    # fr_slindis_l = []
    # fd_slindis_l = []
    # fr_nslindis_l = []
    # fd_nslindis_l = []
    # for i_mc in range(len_mc_l):
    #     if max(sli_l[:,i_mc]) >= cd and max(dis_l[:,i_mc]) >= ca:
    #         fr_slidis_l.append(fr_l[i_mc])
    #         fd_slidis_l.append(fd_l[i_mc])
    #     if max(sli_l[:,i_mc]) < cd and max(dis_l[:,i_mc]) < ca:
    #         fr_nslindis_l.append(fr_l[i_mc])
    #         fd_nslindis_l.append(fd_l[i_mc])
    #     if max(sli_l[:,i_mc]) >= cd and max(dis_l[:,i_mc]) < ca:
    #         fr_slindis_l.append(fr_l[i_mc])
    #         fd_slindis_l.append(fd_l[i_mc])
    # fr_slidis_l = array(fr_slidis_l)
    # fd_slidis_l = array(fd_slidis_l)
    # fr_nslindis_l = array(fr_nslindis_l)
    # fd_nslindis_l = array(fd_nslindis_l)
    #             for i_arr in range(nx):
    #                 sq_dist = (new_arr_x_l[i_res] - arr_x_l[i_arr])**2.0 + \
    #                             (new_arr_y_l[j_res] - arr_y_l[i_arr])**2.0
                    # print(sq_dist)
    #                 if sq_dist <= radius**2.0:
    #                     m+=1
    #             dens[i_res,j_res] = m
        #print(shape(dens))
        #print(dens[0,:])
    #     for ii_res in range(res):
    #         for jj_res in range(fac_res*res):
         #       print(dens[ii_res,jj_res])
    #             if dens[ii_res,jj_res] < 2.:
    #                 print(ii_res,jj_res)
    #                 dens_output[ii_res]=(new_arr_y_l[jj_res])
    #                 break
    # if type == 'min':
    #     for i_res in range(res):
    #         m=0
    #         for j_res in range(fac_res*res):
    #             for i_arr in range(nx):
    #                 sq_dist = (new_arr_x_l[i_res] - arr_x_l[i_arr])**2.0 + \
    #                             (new_arr_y_l[j_res] - arr_y_l[i_arr])**2.0
                    # print(sq_dist)
    #                 if sq_dist <= radius**2.0:
    #                     m+=1
    #             dens[i_res,j_res] = m
        #print(shape(dens))
        #print(dens[0,:])
    #     for ii_res in range(res):
    #         for jj_res in reversed(range(fac_res*res)):
         #       print(dens[ii_res,jj_res])
    #             if dens[ii_res,jj_res] < 1.:
    #                 print(ii_res,jj_res)
    #                 dens_output[ii_res]=(new_arr_y_l[jj_res])
    #                 break
    # print(dens_output)
    # return new_arr_x_l,dens_output
    
def trans_fr(fr_l,fu_l,froude,d_fd):
    n_fdr_l = len(fr_l)
    fu = []
    for i in range(n_fdr_l):
        if fr_l[i] >= froude - d_fd and fr_l[i] <= froude + d_fd:
            fu.append(fu_l[i])
    fu = asarray(fu)
    
    if len(fu)== 0:
        tt = amax(fu_l)
    else:
        tt = amin(fu)     
    return(tt)
                    
                    
def find_boundary(arr_x_l,arr_y_l,res,type='max'):
    nx = len(arr_x_l)
    fac_res =1
    radius = 0.01
    new_arr_x_l = linspace(min(arr_x_l),max(arr_x_l),res)
    new_arr_y_l = linspace(min(arr_y_l),max(arr_y_l),fac_res*res)
    dens = zeros([res,fac_res*res])
    dens_output = zeros(res)
    if type == 'max':
        for i_res in range(res):
            m=0
            for j_res in reversed(range(fac_res*res)):
                for i_arr in range(nx):
                    sq_dist = (new_arr_x_l[i_res] - arr_x_l[i_arr])**2.0 + \
                                (new_arr_y_l[j_res] - arr_y_l[i_arr])**2.0
                    # print(sq_dist)
                    if sq_dist <= radius**2.0:
                        m+=1
                dens[i_res,j_res] = m
        #print(shape(dens))
        #print(dens[0,:])
        for ii_res in range(res):
            for jj_res in range(fac_res*res):
         #       print(dens[ii_res,jj_res])
                if dens[ii_res,jj_res] < 2.:
                    print(ii_res,jj_res)
                    dens_output[ii_res]=(new_arr_y_l[jj_res])
                    break
    if type == 'min':
        for i_res in range(res):
            m=0
            for j_res in range(fac_res*res):
                for i_arr in range(nx):
                    sq_dist = (new_arr_x_l[i_res] - arr_x_l[i_arr])**2.0 + \
                                (new_arr_y_l[j_res] - arr_y_l[i_arr])**2.0
                    # print(sq_dist)
                    if sq_dist <= radius**2.0:
                        m+=1
                dens[i_res,j_res] = m
        #print(shape(dens))
        #print(dens[0,:])
        for ii_res in range(res):
            for jj_res in reversed(range(fac_res*res)):
         #       print(dens[ii_res,jj_res])
                if dens[ii_res,jj_res] < 1.:
                    print(ii_res,jj_res)
                    dens_output[ii_res]=(new_arr_y_l[jj_res])
                    break
    print(dens_output)
    return new_arr_x_l,dens_output

print()
print("\t\t Analyzing data")

def motion_areas():
    data_sli = loadtxt("../output_model_sliding.dat")#,delimiter="")
    data_dis = loadtxt("../output_model_rotation.dat")#,delimiter=",")
    condi = loadtxt("../output_model_conditions.dat")#,delimiter=",")
    nx,ny = shape(data_sli)
    print(nx,ny)
    fd = condi[:,1]
    fr = condi[:,0]
    time = data_sli[:,0]
    n_real = ny-1
    critical_angle = 0.5 * pi
    critical_dista = 3.0
    u_x = []
    u_y = []
    fr_slidis ,fd_slidis,fr_nslindis ,fd_nslindis,fr_slindis ,fd_slindis = \
        motion_space(fr,fd,data_sli[:,1:],data_dis[:,1:],critical_dista,critical_angle)
    
    fr_n, fd_n,dis_n,sli_n = find_boundary_impr(fr,fd,data_dis[:,1:],data_sli[:,1:],50)
    deci1 = zeros([len(fr_n),len(fd_n)])
    print(len(fr_n),len(fd_n))
    for i in range(len(fr_n)):
        for j in range(len(fd_n)):
            print(i,j,fr_n[i],fd_n[j],dis_n[i,j],sli_n[i,j])
            if dis_n[i,j] < 0.5*pi and sli_n[i,j] < 3.0:
                deci1[i,j] = 0.0
            if dis_n[i,j] < 0.5 * pi and sli_n[i,j] >= 3.0:
                deci1[i,j] = 1.0
            if dis_n[i,j] >= 0.5 * pi  and sli_n[i,j] >= 3.0:
                deci1[i,j] = 2.0
    print(deci1)
    plt.pcolor(deci1) 
    plt.show() 


#    print(shape(z_n))
#        find_boundary_impr(fr_nslindis,fd_nslindis,20)
#    b2_x,b2_y = find_boundary_impr(fr_slidis,fd_slidis,20,'min')
    # b1_x = array(b1_x)
    # b1_y = array(b1_y)
    # b2_x = array(b2_x)
    # b2_y = array(b2_y)
    # 
    # p1_x = zeros(len(b1_x))
    # p1_y = zeros(len(b1_x))
    # p1_x[:] = (b1_x[:])
    # p1_y[:] = (b1_y[:])
    # p1_x = append(p1_x,b1_x[-1])
    # p1_y = append(p1_y,0.0)
    # p1_x = append(p1_x,0.0)
    # p1_y = append(p1_y,0.0) 
    # p1_x = append(p1_x,0.0)
    # p1_y = append(p1_y,b1_y[0]) 
    # 
    # p2_x = zeros(len(b1_x))
    # p2_y = zeros(len(b1_x))
    # p2_x[:] = (b1_x[:])
    # p2_y[:] = (b1_y[:])
    # p2_x = append(p2_x,b2_x[-1])
    # p2_y = append(p2_y,b2_y[-1])
    # p2_x = concatenate((p2_x,b2_x[::-1]))
    # p2_y = concatenate((p2_y,b2_y[::-1]))
    # p2_x = append(p2_x,b2_x[0])
    # p2_y = append(p2_y,20.0)
    # p2_x = append(p2_x,b1_x[0])
    # p2_y = append(p2_y,20.0)
    # 
    # p3_x = zeros(len(b2_x))
    # p3_y = zeros(len(b2_x))
    # p3_x[:] = (b2_x[:])
    # p3_y[:] = (b2_y[:])
    # p3_x = append(p3_x,b2_x[-1])
    # p3_y = append(p3_y,20.0)
    # p3_x = append(p3_x,b2_x[0])
    # p3_y = append(p3_y,20.0) 
    # 
    # al = 0.9
    # 
    # fig, ax1 = plt.subplots(nrows=1,ncols=1,figsize=(12,10))
    # font0 = FontProperties()
    # font = font0.copy()
    # font.set_weight('bold')
    # 
    # ax1.tick_params(which='both', width=3)
    # ax1.tick_params(which='major', length=12)
    # ax1.scatter(fr_slidis,fd_slidis,color = 'r',alpha=al)
    # ax1.scatter(fr_nslindis,fd_nslindis,color = 'b',alpha=al)
    # ax1.scatter(fr_slindis,fd_slindis,color = 'g',alpha=al)
    # ax1.fill(p1_x,p1_y,'b',alpha=0.3)
    # ax1.fill(p2_x,p2_y,'g',alpha=0.3)
    # ax1.fill(p3_x,p3_y,'r',alpha=0.3)
    # ax1.set_xlim(min(fr),max(fr))    
    # ax1.set_ylim(min(fd),max(fd))
    # ax1.set_xlabel('Froude Number',fontproperties=font);
    # ax1.set_ylabel('Flow Depth',fontproperties=font);
    # plt.show()
    # 
    # make_outputfile(fig,tt='test',res=200)

def simple_plot():
    aaa,bbb,ccc,MC_set,n_proc,rhos_min,rhos_max,fd_min,fd_max,\
        fr_min,fr_max,cd_min,cd_max,cl_min,cl_max,mu_min,mu_max,\
        cma_min,cma_max,slope,delta,num_time,end_time=read_input('../boulder.in')
    if delta == 0.0:
        data_sli = loadtxt("../output_model_sliding.dat")#,delimiter=",")
        data_dis = loadtxt("../output_model_rotation.dat")#,delimiter=",")
        condi = loadtxt("../output_model_conditions.dat")#,delimiter=",")
        nx,ny = shape(data_sli)
        print(nx,ny)
        fd = condi[:,1]
        fr = condi[:,0]
        time = data_sli[:,0]
        n_real = ny-1
        critical_angle = 0.5 * pi
        critical_dista = 3.0
        u_x = []
        u_y = []
        fr_slidis ,fd_slidis,fr_nslindis ,fd_nslindis,fr_slindis ,fd_slindis = \
            motion_space_c(fr,fd,data_sli[:,1:],data_dis[:,1:],critical_dista,critical_angle)
        fig, ax1 = plt.subplots(nrows=1,ncols=1,figsize=(12,10))
        font0 = FontProperties()
        font = font0.copy()
        font.set_weight('bold')
        al = 1.0
        ax1.tick_params(which='both', width=3)
        ax1.tick_params(which='major', length=12)
        ax1.scatter(fr_slidis,fd_slidis,color = 'r',alpha=al)
        ax1.scatter(fr_nslindis,fd_nslindis,color = 'b',alpha=al)
        ax1.scatter(fr_slindis,fd_slindis,color = 'g',alpha=al)
        ax1.set_xlim(min(fr),max(fr))    
        ax1.set_ylim(min(fd),max(fd))
        ax1.set_xlabel('Froude Number',fontproperties=font);
        ax1.set_ylabel('Flow Depth',fontproperties=font);
        plt.show()
    if delta > 0.0:
        data_dis = loadtxt("../output_model_rotation.dat")#,delimiter=",")
        condi = loadtxt("../output_model_conditions.dat")#,delimiter=",")
        nx,ny = shape(data_dis)
        print(nx,ny)
        fd = condi[:,1]
        fr = condi[:,0]
        time = data_dis[:,0]
        n_real = ny-1
        critical_angle = 0.5 * pi
        critical_dista = 3.0
        u_x = []
        u_y = []
        fr_dis ,fd_dis,fr_ndis ,fd_ndis = \
            motion_space_p(fr,fd,data_dis[:,1:],critical_angle)
        fig, ax1 = plt.subplots(nrows=1,ncols=1,figsize=(12,10))
        font0 = FontProperties()
        font = font0.copy()
        font.set_weight('bold')
        al = 1.0
        ax1.tick_params(which='both', width=3)
        ax1.tick_params(which='major', length=12)
        ax1.scatter(fr_dis,fd_dis,color = 'r',alpha=al)
        ax1.scatter(fr_ndis,fd_ndis,color = 'b',alpha=al)
        ax1.set_xlim(min(fr),max(fr))    
        ax1.set_ylim(min(fd),max(fd))
        ax1.set_xlabel('Froude Number',fontproperties=font);
        ax1.set_ylabel('Flow Depth',fontproperties=font);
        plt.show()
    make_outputfile(fig,tt='simple',res=200)
 
    
def para_simple():
    sl_change = [0.0]
#    dl_change = [0.1, 0.5,1.0,1.5, 2.0, 2.5,3.0, 3.5,4.0]
    dl_change = [0.1]+ list(linspace(0.5,4.0,15))
    m=0
    al = 1.0
    trans = []
    froude_value = 1.0
    df = 0.02
    plt.figure(figsize=(20,20))
    for i in range(len(sl_change)):
        for j in range(len(dl_change)):
            m = m +1
            filename = 'out_'+det_runname(sl_change[i],dl_change[j])
            fl1 = '../'+filename+'_conditions.dat'
            fl2 = '../'+filename+'_rotation.dat'
            data_dis = loadtxt(fl2)#,delimiter=",")
            condi = loadtxt(fl1)#,delimiter=",")
            fd = condi[:,1]
            fr = condi[:,0]
            time = data_dis[:,0]
            critical_angle = 0.5 * pi
            critical_dista = 3.0
            fr_dis ,fd_dis,fr_ndis ,fd_ndis = \
            motion_space_p(fr,fd,data_dis[:,1:],critical_angle)
            # plt.subplot(3,3,m)
            # plt.scatter(fr_dis,fd_dis,color = 'r',alpha=al)
            # plt.scatter(fr_ndis,fd_ndis,color = 'b',alpha=al)
            trans.append(trans_fr(fr_dis,fd_dis,froude_value,df))
            # plt.plot([froude_value],[trans[-1]],'go')
            # plt.plot([froude_value-df,froude_value+df],[trans[-1],trans[-1]],'g-',lw=10)
#            plt.axvline([0.6])
#            print(m,trans[-1])
    plt.scatter(dl_change,trans,color='k')
    plt.show()
    

if __name__ == '__main__':
#    para_simple()
        
    simple_plot()
#    motion_areas()
