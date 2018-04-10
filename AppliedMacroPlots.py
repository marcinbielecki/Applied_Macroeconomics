# Marcin Bielecki, University of Warsaw, mbielecki@wne.uw.edu.pl

import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rcParams
rcParams['axes.autolimit_mode'] = 'round_numbers'
rcParams['axes.xmargin'] = 0
rcParams['axes.ymargin'] = 0

rcParams['legend.frameon'] = False

def NiceFigs(fig, ax, northeast=True):
    # below is algorithm for nice axes
 
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
 
    # removing the default axis on all sides:
    for side in ['bottom','right','top','left']:
        ax.spines[side].set_visible(False)
 
    # removing the axis ticks
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
 
    # get width and height of axes object to compute
    # matching arrowhead length and width
    dps = fig.dpi_scale_trans.inverted()
    bbox = ax.get_window_extent().transformed(dps)
    width, height = bbox.width, bbox.height
     
    # manual arrowhead width and length
    hw = 1./50. * (ymax-ymin)
    hl = 1./50. * (xmax-xmin)
    ohg = 0.3 # arrow overhang
     
    # compute matching arrowhead length and width
    yhw = hw/(ymax-ymin) * (xmax-xmin) * height/width
    yhl = hl/(xmax-xmin) * (ymax-ymin) * width/height
    
    # axis line width
    lw = 1./1000. * (ymax-ymin)
    ylw = lw/(ymax-ymin) * (xmax-xmin) * height/width
     
    # draw x and y axis
    if northeast==True:
        ax.arrow(xmin, ymin, xmax-xmin, 0., fc='k', ec='k', width=lw,
                 head_width=hw, head_length=hl, overhang=ohg,
                 length_includes_head=True, clip_on=False, zorder=100)
         
        ax.arrow(xmin, ymin, 0., ymax-ymin, fc='k', ec='k', width=ylw,
                 head_width=yhw, head_length=yhl, overhang=ohg,
                 length_includes_head=True, clip_on=False, zorder=100)
    else:
        ax.arrow(xmin, 0, xmax-xmin, 0., fc='k', ec='k', width=lw,
                 head_width=hw, head_length=hl, overhang=ohg,
                 length_includes_head=True, clip_on=False, zorder=100)
         
        ax.arrow(0, ymin, 0., ymax-ymin, fc='k', ec='k', width=ylw,
                 head_width=yhw, head_length=yhl, overhang=ohg,
                 length_includes_head=True, clip_on=False, zorder=100)
     
    # clip_on = False if only positive x or y values.
    
    return None

    
def PlotBudgetSet(y_0=1, y_1=1, r=0.05, figsize=(5, 5), fig_ax=None, full_labels=False):

    c_0_max = y_0 + y_1/(1+r)
    c_1_max = (1+r)*c_0_max
    c_0_mid = (y_0 + y_1/(1+r))/2
    c_1_mid = (1+r)*(y_0-c_0_mid) + y_1

    c_0 = np.linspace(0, c_0_max)

    c_max = np.ceil(max(c_0_max, c_1_max))

    if fig_ax == None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = fig_ax[0]
        ax = fig_ax[1]

    p = ax.plot(c_0, (1+r)*(y_0-c_0) + y_1, linewidth=2, label='Budget line')
    color = p[0].get_color()

    if full_labels == True:
        ax.fill_between(c_0, 0, (1+r)*(y_0-c_0) + y_1, alpha=0.25, label='Budget set')
    else:
        ax.fill_between(c_0, 0, (1+r)*(y_0-c_0) + y_1, alpha=0.25)
    ax.plot(y_0, y_1, 'o', color=color)

    ax.vlines(y_0, 0, y_1, color=color, linestyle='dashed', linewidth=1)
    ax.hlines(y_1, 0, y_0, color=color, linestyle='dashed', linewidth=1)

    ax.set_xlim(0, c_max)
    ax.set_ylim(0, c_max)

    NiceFigs(fig, ax)

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    if full_labels==True:
        plt.xticks((0, y_0, c_0_max, xmax), 
                   ('0', '$y_{t}$', '$y_{t}+y_{t+1}/(1+r)$', '$c_{t}$'), 
                   fontsize=12)
        plt.yticks((0, y_1, c_1_max, ymax), 
                   ('', '$y_{t+1}$', '$(1+r)\cdot y_{t}+y_{t+1}$', '$c_{t+1}$'), 
                   fontsize=12)
    
        ax.annotate('Slope = $-(1+r)$', va='top', ha='right', 
                    xy=(c_0_mid, c_1_mid), xytext=(c_0_max, c_1_max),
                    arrowprops=dict(arrowstyle='->', 
                                    connectionstyle='arc3, rad=-0.2'), 
                    fontsize=12)
        
        l = plt.legend(loc='upper right', fontsize=12, frameon=True)
        l.get_frame().set_linewidth(0)
        l.get_frame().set_alpha(1)
        
    else:
        plt.xticks((0, y_0, xmax), ('0', '$y_{t}$', '$c_{t}$'), fontsize=12)
        plt.yticks((0, y_1, ymax), ('', '$y_{t+1}$', '$c_{t+1}$'), fontsize=12)

    # plt.savefig('test.pdf')
    return None


def PlotIndifferenceCurves(y_0=1, y_1=1, β=0.95, r=0.05, levels=[1], figsize=(5, 5), fig_ax=None, full_labels=False):

    c_0_max = y_0 + y_1/(1+r)
    c_1_max = (1+r)*c_0_max

    c_max = np.ceil(max(c_0_max, c_1_max))

    c_0 = np.linspace(0.01, c_max, 1e3)
    
    if fig_ax == None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = fig_ax[0]
        ax = fig_ax[1]

    for level in levels:
        # ax.plot(c_0, (level/c_0)**(1/β), linewidth=2, label='U = '+str(level))
        if full_labels == True:
            ax.plot(c_0, np.exp((level-np.log(c_0))/β), linewidth=2, label='U = {:.1f}'.format(level))
        else:
            ax.plot(c_0, np.exp((level-np.log(c_0))/β), linewidth=2, label='Indifference curve')

    ax.set_xlim(0, c_max)
    ax.set_ylim(0, c_max)
    
    NiceFigs(fig, ax)

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    
    if full_labels == True:
        plt.xticks((0, xmax), ('0', '$c_{t}$'), fontsize=12)
        plt.yticks((0, ymax), ('', '$c_{t+1}$'), fontsize=12)
        l = plt.legend(loc='upper right', fontsize=12, frameon=True)
        l.get_frame().set_linewidth(0)
        l.get_frame().set_alpha(1)

    # plt.savefig('test.pdf')
    return None


def SolveUtilityMaximizationProblem(y_0=1, y_1=1, β=0.95, r=0.05, figsize=(5, 5)):

    c_0_max = y_0 + y_1/(1+r)
    c_0 = np.linspace(0, c_0_max)
    
    c_0_opt = 1/(1+β) * (y_0 + y_1/(1+r))
    c_1_opt = β*(1+r)*c_0_opt
    U_opt = np.log(c_0_opt) + β*np.log(c_1_opt)
    
    fig, ax = plt.subplots(figsize=figsize)

    PlotBudgetSet(y_0, y_1, r, fig_ax=(fig, ax))
    PlotIndifferenceCurves(y_0, y_1, β, levels=[U_opt], fig_ax=(fig, ax))

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    
    p = ax.plot(c_0, β*(1+r)*c_0, linewidth=2, label='Euler equation')
    color = p[0].get_color()
    
    ax.plot(c_0_opt, c_1_opt, 'o', color=color)

    plt.vlines(c_0_opt, 0, c_1_opt, color=color, linestyle='dashed', linewidth=1)
    plt.hlines(c_1_opt, 0, c_0_opt, color=color, linestyle='dashed', linewidth=1)
    
    l = plt.legend(loc='upper right', fontsize=12)
    l.get_frame().set_linewidth(0)
    l.get_frame().set_alpha(1)
    
    if c_0_opt == y_0:
        plt.xticks((0, y_0, xmax), ('0', '$c_{t} = y_{t}$', '$c_{t}$'), fontsize=12)
        plt.yticks((0, y_1, ymax), ('', '$c_{t+1} = y_{t+1}$', '$c_{t+1}$'), fontsize=12)
    elif abs(c_0_opt-y_0) < xmax/20:
        if c_0_opt > y_0:
            plt.xticks((0, (y_0+c_0_opt)/2, xmax), ('0', '$c_{t} > y_{t}$', '$c_{t}$'), fontsize=12)
            plt.yticks((0, (y_1+c_1_opt)/2, ymax), ('', '$c_{t+1} < y_{t+1}$', '$c_{t+1}$'), fontsize=12)
        else:
            plt.xticks((0, (y_0+c_0_opt)/2, xmax), ('0', '$c_{t} < y_{t}$', '$c_{t}$'), fontsize=12)
            plt.yticks((0, (y_1+c_1_opt)/2, ymax), ('', '$c_{t+1} > y_{t+1}$', '$c_{t+1}$'), fontsize=12)
    else:
        plt.xticks((0, y_0, c_0_opt, xmax), ('0', '$y_{t}$', '$c_{t}$', '$c_{t}$'), fontsize=12)
        plt.yticks((0, y_1, c_1_opt, ymax), ('', '$y_{t+1}$', '$c_{t+1}$', '$c_{t+1}$'), fontsize=12)

    # plt.savefig('test.pdf')
    
    print("Consumer's income:")
    print('  y_0: {:.4f}'.format(y_0))
    print('  y_1: {:.4f}'.format(y_1))
    print('')
    
    print('Optimal consumption:')
    print('  c_0: {:.4f}'.format(c_0_opt))
    print('  c_1: {:.4f}'.format(c_1_opt))
    print('')
    
    if c_0_opt > y_0:
        print('c_0 > y_0: Consumer is a borrower')
    elif c_0_opt < y_0:
        print('c_0 < y_0: Consumer is a saver')
    else:
        print('c_0 = y_0: Consumer is neither a borrower nor a saver')
    print('')
        
    return None


def ShowIncomeAndSubstitutionEffects(y_0=1, y_1=1, β=0.95, r=0.05, r_new=1, figsize=(5, 5)):
    c_0_opt = 1/(1+β) * (y_0 + y_1/(1+r))
    c_1_opt = β*(1+r)*c_0_opt
    U_opt = np.log(c_0_opt) + β*np.log(c_1_opt)
    
    c_0_new = 1/(1+β) * (y_0 + y_1/(1+r_new))
    c_1_new = β*(1+r_new)*c_0_new
    U_new = np.log(c_0_new) + β*np.log(c_1_new)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    c_0_max = y_0 + y_1/(1+r)
    c_0_max_new = y_0 + y_1/(1+r_new)
    
    c_1_max = (1+r)*c_0_max
    c_1_max_new = (1+r_new)*c_0_max_new

    c_max = np.ceil(max(c_0_max, c_0_max_new, c_1_max, c_1_max_new))
    
    c_0 = np.linspace(0, c_max)
    
    ax.plot(c_0, (1+r)*(y_0-c_0) + y_1, linewidth=1, color='C0')
    ax.plot(c_0, (1+r_new)*(y_0-c_0) + y_1, linewidth=2, color='C0')
    
    ax.plot(c_0[1:], np.exp((U_opt-np.log(c_0[1:]))/β), linewidth=1, color='C1')
    ax.plot(c_0[1:], np.exp((U_new-np.log(c_0[1:]))/β), linewidth=2, color='C1')
    
    c_0_sub = np.exp((U_opt-β*np.log(β*(1+r_new)))/(1+β))
    c_1_sub = β*(1+r_new)*c_0_sub
    
    ax.plot(c_0, (1+r_new)*(c_0_sub + c_1_sub/(1+r_new) - c_0), linewidth=1, color='k')

    ax.set_xlim(0, c_max)
    ax.set_ylim(0, c_max)

    NiceFigs(fig, ax)

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    
    ax.plot(y_0, y_1, 'o', color='C0', label='Initial endowment, $E$')
    
    ax.plot(c_0_opt, c_1_opt, 'o', color='C2', markerfacecolor='none', label='Optimal choice for old $r$, $O$')

    plt.vlines(c_0_opt, 0, c_1_opt, color='C2', linestyle='dashed', linewidth=1)
    plt.hlines(c_1_opt, 0, c_0_opt, color='C2', linestyle='dashed', linewidth=1)
    
    ax.plot(c_0_new, c_1_new, 'o', color='C2', label='Optimal choice for new $r$, $O^{\prime}$')

    plt.vlines(c_0_new, 0, c_1_new, color='C2', linestyle='dashed', linewidth=1)
    plt.hlines(c_1_new, 0, c_0_new, color='C2', linestyle='dashed', linewidth=1)
    
    ax.plot(c_0_sub, c_1_sub, 'o', color='k', label='Pure substitution effect, $S$')
    
    l = plt.legend(loc='upper right', fontsize=12)
    l.get_frame().set_linewidth(0)
    l.get_frame().set_alpha(1)
    
    if c_0_opt == c_0_new:
        plt.xticks((0, c_0_opt, xmax), ('0', '$c_{t} = c^{\prime}_{t}$', '$c_{t}$'), fontsize=12)
    elif abs(c_0_opt-c_0_new) < xmax/20:
        if c_0_opt > c_0_new:
            plt.xticks((0, (c_0_new+c_0_opt)/2, xmax), ('0', '$c^{\prime}_{t} < c_{t}$', '$c_{t}$'), fontsize=12)
        else:
            plt.xticks((0, (c_0_new+c_0_opt)/2, xmax), ('0', '$c^{\prime}_{t} > c_{t}$', '$c_{t}$'), fontsize=12)
    else:
        plt.xticks((0, c_0_new, c_0_opt, xmax), ('0', '$c^{\prime}_{t}$', '$c_{t}$', '$c_{t}$'), fontsize=12)
    
    if c_1_opt == c_1_new:
        plt.yticks((0, c_1_opt, ymax), ('', '$c_{t+1} = c^{\prime}_{t+1}$', '$c_{t+1}$'), fontsize=12)
    elif abs(c_1_opt-c_1_new) < ymax/20:
        if c_1_opt > c_1_new:
            plt.yticks((0, (c_1_new+c_1_opt)/2, ymax), ('', '$c^{\prime}_{t+1} > c_{t+1}$', '$c_{t+1}$'), fontsize=12)
        else:
            plt.yticks((0, (c_1_new+c_1_opt)/2, ymax), ('', '$c^{\prime}_{t+1} < c_{t+1}$', '$c_{t+1}$'), fontsize=12)
    else:
        plt.yticks((0, c_1_new, c_1_opt, ymax), ('', '$c^{\prime}_{t+1}$', '$c_{t+1}$', '$c_{t+1}$'), fontsize=12)
    
    return None


def SolveUtilityMaximizationProblemCRRA(y_0=1, y_1=1, β=1, σ=1, r=0, fig_ax=None):
    C_1_max = y_0 + y_1/(1+0)
    C_2_max = (1+1)*C_1_max

    C_max = np.ceil(max(C_1_max, C_2_max))/3
    
    PDV_Y = y_0 + y_1/(1+r)
    
    C_1_opt = 1/(1 + β**(1/σ) * (1+r)**(1/σ-1)) * PDV_Y
    C_2_opt = (β * (1+r))**(1/σ) * C_1_opt
    
    fig = fig_ax[0]
    ax = fig_ax[1]
    
    ax.plot(C_1_opt, C_2_opt, 'o', label='r =  {:.1f}'.format(r))
    
#     ax.set_xlim(0, C_max)
#     ax.set_ylim(0, C_max)
    
#     NiceFigs(fig, ax)

#     xmin, xmax = ax.get_xlim()
#     ymin, ymax = ax.get_ylim()
    
#     plt.xticks((0, xmax), ('', '$C_{1}$'), fontsize=12)
#     plt.yticks((0, ymax), ('', '$C_{2}$'), fontsize=12)
    
    plt.legend(frameon=False)
    
#     print('Optimal consumption:')
#     print('  C_1: {:.4f}'.format(C_1_opt))
#     print('  C_2: {:.4f}'.format(C_2_opt))
#     print('')
    
    return None
