import scipy.io as sio
import numpy as np
import numpy.ma as ma
import seaborn as sns
import pandas as pd
from pathlib import Path
import itertools
import matplotlib.pyplot as plt
from PIL import Image
import matplotlib.image as mpimg
from scipy.ndimage import gaussian_filter
from scipy import stats
from mpl_toolkits.axes_grid1 import make_axes_locatable
from skimage.transform import resize
import seaborn as sns
import os
import re
import cv2
from scipy.stats.stats import pearsonr
from matplotlib.colors import LogNorm
import pickle

path = "/Users/terlau/eye_data/"
pic_path = "/Users/terlau/stimuli_640x480/"
shown_imgages = np.loadtxt('#_of_pics.dat',delimiter=',')
#img_33 = mpimg.imread(pic_path)

def get_viewed_img(pic_path):
    stimulus = []
    for subdir, dirs, files in os.walk(pic_path):
        for file in files:
            if file == '.DS_Store' or file.endswith("_35pxl.png") or file.endswith("_35pxl.bmp"):
                continue
            if file.endswith(".bmp") or file.endswith(".png"):
                number = re.findall(r'\d+', file)
                file_number = int(number[0])

                if file_number in shown_imgages:
                    image = os.path.join(subdir, file) #Image.open
                    #print(image)
                    stimulus.append(image)
    return stimulus

stimuli_paths = get_viewed_img(pic_path)
print(type(stimuli_paths[0]), len(stimuli_paths), stimuli_paths[0])

def get_data(data_path):

    pic_number = []
    category = []
    subjects = []
    x_coord = []
    y_coord = []
    images = []
    blinks = []
    saccades = []
    deep_gaze_files = []
    for subdir, dirs, files in os.walk(path):
        for file in files:
            if file == '.DS_Store' or file == 'eyedata.mat':
                continue
            if file.endswith(".mat"):
                #print("File:", file)
                sub = re.findall(r'\d+', file)[0] #'^[0-9]+'  search gives only re.match object which is not iterable
                #print("SUBJECT", int(sub))
                data = sio.loadmat(os.path.join(subdir, file))
                # this corresponds to trial in the .mat file, shape (150,6,1301)
                sub_trial = data['timelock'][0][0][3]
                #print("Sub Trial Channels 4 and 5", sub_trial.shape, sub_trial[:,[4,5]])

                sub_trialinfo = data['timelock'][0][0][4]
                #time = data['timelock'][0][0][0]
                #print(type(sub_trial), sub_trial.shape)
                if len(sub_trial[:,0]) != 150:
                    #print(file) #eye_sub-63.mat, eye_sub-60.mat, eye_sub-10.mat are excluded due to missing data for images
                    continue
                else:
                    # is second and third column correct vor trial (viewing behavior) as well?!
                    subTrialinfo = sub_trialinfo[:,[1,2]]
                    #print(sub_trial[:, [1, 2], :].shape)
                    #trial_images.append(sub_trial[:, [1, 2], :])
                    for i in range(0,150):
                        subjects.append(sub)
                        images.append(i)
                        category.append(subTrialinfo[i][0])
                        pic_number.append(subTrialinfo[i][1])
                        img = sub_trial[i]
                        img_i = img[[1, 2, 4, 5], :]
                        img_fin = img_i[:,300:801]
                        #print(img_fin.shape)
                        blinks.append(img_fin[2])
                        saccades.append(img_fin[3])
                        #plt.plot(img_fin[0], img_fin[1])
                        #plt.show()

                        x_coord.append(img_fin[0])
                        y_coord.append(img_fin[1])


                    # dictionary of lists
                    dicto = {'Subjects': subjects, 'Images': images, 'Category': category, 'pic_#': pic_number, 'x': x_coord, 'y': y_coord, 'blinks': blinks, 'sacc': saccades}
                    df = pd.DataFrame(dicto)

            else:
                #print(file)
                deep_gaze_files.append(file)
    return df, deep_gaze_files



    #return trial_images, trialinfo, subjects

data, deep_gaze = get_data(path)
data["Subjects"] = data['Subjects'].astype('int')


#print("DTypes", data.dtypes)
print("Data", data)
print("Deep Gaze", len(deep_gaze))
#print(data.loc[0:150]['Subjects'])
# Length of columns = 13200 (88*150)
print(data.loc[:]['sacc'].shape)

def get_img_path(name):
    cat = {
    "fractals" : 1,
    "landscapes" : 2,
    "naturals1" : 3,
    "streets1" : 4,
    "streets2" : 5
    }
    #have to change tuple to list to allow for value assignment later on
    lst = list(name)

    #print("Categories", cat)
    for key, value in cat.items():
        if lst[0] == value:
            # if category of pic is matched in dictonary assign key (cat name) to first value
            lst[0]=key

    pic =  int(lst[1])
    if lst[0] in ('landscapes', 'streets2'):
        path_img = os.path.join(pic_path, lst[0], 'image'+str(pic)+'.png')
    else:
        path_img = os.path.join(pic_path, lst[0], str(pic)+'.bmp')

    return path_img, lst[0]


cat = {
"fractals" : 1,
"landscapes" : 2,
"naturals1" : 3,
"streets1" : 4,
"streets2" : 5
}

data_image = data.groupby(['Category', 'pic_#'])#, as_index=False)
#print("Data Grouped per image:", data_image.head())
#data_sub = data.groupby('Subjects')
#print("Data Grouped per subject:", data_sub.head())

# 3x3 renders a 2x2, so gotta use 4x4
x_edges = np.linspace(192, 832, num=4) # 455+1 since shape of hist otherwise 454*340
y_edges = np.linspace(144, 624, num=4) #342


n = 0
corr_values = []
rho_values = []
dataframes_sub = []
pearson_mean = []
rho_mean = []


for name, group in data_image:
    n += 1
    # group is per subject, so 150 rows, 8 columns

    print("Type", group.head(), group, type(group))
    #print("Current Image", group['Category'], group['pic_#']) --> 150*same category, 150*pic numbers (these are obviously different then)

    x_coord_img = group["x"].mean()
    y_coord_img = group["y"].mean()
    blinks = group['blinks'].mean()
    sacc = group['sacc'].mean()

    #print("Coords averaged per image:", x_coord_img.shape, y_coord_img.shape)
    #print("Blink and saccades averaged per image:", blinks.shape, sacc.shape)

    #print(sacc)

    sacc_bi = (sacc > 0.5).astype(int)
    blinks_bi = (blinks > 0.5).astype(int)
    #print("Unique values blinks:", np.unique(blinks_bi))
    #print("Unique values saccades:", np.unique(sacc_bi))

    mask = (sacc_bi == 0) & (blinks_bi == 0)
    print("Mask:", mask, np.unique(mask))
    x_coord_img_t = x_coord_img[mask]
    y_coord_img_t = y_coord_img[mask]

    #marker_size=7
    #plt.scatter(x_coord_img, y_coord_img,  marker_size, alpha=0.8, label='viewing')
    #plt.scatter(x_coord_img_t, y_coord_img_t, marker_size, alpha= 0.5,label='masked viewing')
    #plt.legend()
    #plt.show()
    x_edg = np.linspace(192, 832, num=455)
    y_edg = np.linspace(144, 624, num=342)
    h_all_bins, xedg, yedg = np.histogram2d(x_coord_img_t, y_coord_img_t,  bins=[x_edg, y_edg])
    h_all = h_all_bins.T
    #H_flipped = np.flipud(H)
    H_smoothed = gaussian_filter(h_all, sigma=10)


    hist, xedges, yedges = np.histogram2d(x_coord_img_t, y_coord_img_t,  bins=[x_edges, y_edges])
    #print("Hist", type(hist), hist.shape, hist)
    H = hist.T
    #H_flipped = np.flipud(H)
    #print("Histogram flipped:", H)
    #H_smoothed_5 = gaussian_filter(H, sigma=5)
    #H_smoothed_10 = gaussian_filter(H, sigma=10)

    #get image that was shown for the averaged viewing behavior
    img_path, pic_cat = get_img_path(name)
    print(img_path)
    #save them
    #img_5 = Image.fromarray(H_smoothed_5)
    #img_5=Image.fromarray((H_smoothed_5 * 255).astype('uint8'), mode='L')
    #img_5.save("FDMs/" + str(pic_cat) + '_' + str(name[1]) + 'sigma_5' + '.png')

    #fdm_5= np.save("FDMs/" + str(pic_cat) + '_' + str(name[1]) + 'sigma_5', H_smoothed_5)    # .npy extension is added if not given

    #img_10=Image.fromarray((H_smoothed_10 * 255).astype('uint8'), mode='L')
    #img_10.save("FDMs/" + str(pic_cat) + '_' + str(name[1]) + 'sigma_10' + '.png')

    #fdm_10= np.save("FDMs/" + str(pic_cat) + '_' + str(name[1]) + 'sigma_10', H_smoothed_10)


    cat = os.path.split(os.path.split(img_path)[0])[1]
    sub_list = [x for x in deep_gaze if cat in x]
    #print(str(int(name[1])))
    deep_gaze_current_file = [x for x in sub_list if "_"+str(int(name[1]))+"." in x]
    print(deep_gaze_current_file[0], img_path)
    deep_gaze_path = os.path.join(path, deep_gaze_current_file[0])
    img = mpimg.imread(img_path)
    #deep_gaze_II = mpimg.imread(deep_gaze_path)
    deep_gaze_II = np.load(deep_gaze_path)
    #flip deep gaze array
    #deep_gaze_II = D.T
    #print("Deep Gaze", type(deep_gaze_II), deep_gaze_II.shape, deep_gaze_II)
    res = cv2.resize(deep_gaze_II, dsize=(3, 3), interpolation=cv2.INTER_AREA) # interpolation=cv2.INTER_AREA recommended for downsampling cv2.INTER_CUBIC
    #print("Rescaled using opencv resize, based on interpolation", res.shape, res)
    res_deep = resize(deep_gaze_II, (3, 3))
    #print("Rescaled using scikit-image resize", res_deep.shape, res_deep) #this might be the right one!
    #for displaying images I need to use nan, for correlation delete!
    res_deep[1,1] = np.nan
    #res[1,1] = np.nan
    #print(res_deep)
    #H_smoothed[1,1] = np.nan
    H[1,1] = np.nan
    #print("Histogram not smoothed", H)
    deep_8 = np.delete(res_deep, 4, None)
    print("Deep Gaze", deep_8, len(deep_8), type(deep_8))
    fdm = np.delete(H, 4, None)
    print("FDM", fdm, len(fdm), type(fdm))

    fig, ax = plt.subplots(2,2)
    ax[0][0].imshow(img, interpolation='nearest', extent=[xedg[0], xedg[-1],yedg[0],yedg[-1]]) #interpolation='gaussian'
    a0 = ax[0][0].imshow(H_smoothed, interpolation='nearest', origin = 'lower', extent=[xedg[0], xedg[-1],yedg[0],yedg[-1]], alpha=0.7) #cmap=plt.cm.RdBu
    ax[0][0].axis('off')
    divider1 = make_axes_locatable(ax[0][0])
    cax1 = divider1.append_axes("right", size="3%", pad=0.05)
    cbar1 = fig.colorbar(a0, ax=ax[0][0], cax = cax1)

    ax[0][1].imshow(img, interpolation='nearest', extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]]) #interpolation='gaussian'
    a1 = ax[0][1].imshow(H, interpolation='nearest', origin = 'lower', extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]], alpha=0.7) #cmap=plt.cm.RdBu
    ax[0][1].axis('off')
    divider2 = make_axes_locatable(ax[0][1])
    cax2 = divider2.append_axes("right", size="3%", pad=0.05)
    cbar2 = fig.colorbar(a1, ax=ax[0][1], cax = cax2)

    ax[1][0].imshow(img, interpolation='nearest', extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]]) #interpolation='gaussian'
    a2 = ax[1][0].imshow(deep_gaze_II, interpolation='nearest', origin = 'lower', extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]], alpha=0.7) #cmap=plt.cm.RdBu
    ax[1][0].axis('off')
    divider3 = make_axes_locatable(ax[1][0])
    cax3 = divider3.append_axes("right", size="3%", pad=0.05)
    cbar3 = fig.colorbar(a2, ax=ax[1][0], cax = cax3)

    ax[1][1].imshow(img, interpolation='nearest', extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]]) #interpolation='gaussian'
    a3 = ax[1][1].imshow(res_deep, interpolation='nearest', extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]], alpha=0.7) #cmap=plt.cm.RdBu
    ax[1][1].axis('off')
    divider4 = make_axes_locatable(ax[1][1])
    cax4 = divider4.append_axes("right", size="3%", pad=0.05)
    cbar4 = fig.colorbar(a3, ax=ax[1][1], cax = cax4)
    plt.show()
    """
    fig, axes = plt.subplots(nrows=2, ncols=3)
    axes[0][0].imshow(img, interpolation='nearest', extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]]) #interpolation='gaussian'
    axes[0][0].imshow(H_smoothed_5, interpolation='nearest', origin = 'lower', extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]], alpha=0.7, vmin=0, vmax=1) #cmap=plt.cm.RdBu
    axes[0][0].set_title('FDM sigma-5 for image {} {}'.format(pic_cat, name[1]))
    axes[0][0].axis('off')

    axes[0][1].imshow(H_smoothed_5, interpolation='nearest', origin = 'lower', extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]], alpha=0.7, vmin=0, vmax=1)
    axes[0][1].set_title('FDM sigma-5 for image {} {}'.format(pic_cat, name[1]))
    axes[0][1].axis('off')

    axes[0][2].imshow(img, interpolation='nearest', extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]]) #interpolation='gaussian'
    axes[0][2].imshow(deep_gaze_II, interpolation='nearest', extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]], alpha=0.7) # cmap=plt.cm.RdBu  origin = 'lower',
    axes[0][2].set_title('DeepGaze prediction for image {} {}'.format(pic_cat, name[1]))
    axes[0][2].axis('off')
    #axes[1].set_xlabel('# of bins x axis')
    #axes[1].set_ylabel('# of bins y axis')
    # scaling of colorbar - 0-1.5 the best?


    axes[1][0].imshow(img, interpolation='nearest', extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]])
    axes[1][0].imshow(H_smoothed_10, interpolation='nearest', origin = 'lower', extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]], alpha=0.7, vmin=0, vmax=1)
    #divider1 = make_axes_locatable(axes[0][1])
    #cax1 = divider1.append_axes("right", size="3%", pad=0.05)
    #cbar1 = fig.colorbar(a0, ax=axes[0][1], cax = cax1)
    axes[1][0].set_title('FDM sigma-10 for image {} {}'.format(pic_cat, name[1]))
    axes[1][0].axis('off')


    axes[1][1].imshow(H_smoothed_10, interpolation='nearest', origin = 'lower', extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]], alpha=0.7, vmin=0, vmax=1)
    #a1=
    #divider2 = make_axes_locatable(axes[1][1])
    #cax2 = divider2.append_axes("right", size="3%", pad=0.05)
    #cbar2 = fig.colorbar(a1, ax=axes[1][1], cax = cax2)
    axes[1][1].set_title('FDM sigma-10 for image {} {}'.format(pic_cat, name[1]))
    axes[1][1].axis('off')



    axes[1][2].imshow(img, interpolation='nearest', extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]]) #interpolation='gaussian'
    axes[1][2].imshow(deep_gaze_II, interpolation='nearest', extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]], alpha=0.7) # cmap=plt.cm.RdBu  origin = 'lower',
    #a1=
    #divider2 = make_axes_locatable(axes[1][1])
    #cax2 = divider2.append_axes("right", size="3%", pad=0.05)
    #cbar2 = fig.colorbar(a1, ax=axes[1][1], cax = cax2)
    axes[1][2].set_title('DeepGaze prediction for image {} {}'.format(pic_cat, name[1]))
    axes[1][2].axis('off')

    #axes[1][2].imshow(img, interpolation='nearest') #interpolation='gaussian'
    #axes[1][2].set_title('Original image {} {}'.format(pic_cat, name[1]))
    #axes[1][2].axis('off')
    #cbar2.ax.set_ylabel('Test', rotation = 270)

    plt.tight_layout()
    #fig.colorbar(test)
    #fig.colorbar(test, ax=axes.ravel().tolist(), shrink=0.5)

    plt.show()
    """
    """
    appended_data = []
    #this loop is only needed if I want values for single subject per image or for single image per subject
    for index, sub_data in group.iterrows():
        #sub_data is a series, it has 5 columns, so it basically gives me the values for each image/subject, one image at a time
        print("SUBDATA", sub_data, type(sub_data))

        x_coord_sub = sub_data["x"]
        y_coord_sub = sub_data["y"]
        #print(x_coord_sub.shape)
        #print(y_coord_sub.shape)

        blinks = sub_data['blinks']
        sacc = sub_data['sacc'] #type numpy ndarray
        #np.where(sacc <= 0.5, 0, 1)
        sacc_bi = (sacc > 0.5).astype(int)
        blinks_bi = (blinks > 0.5).astype(int)
        #print(np.unique(blinks_bi))

        #mask_test = (sacc_bi == 0)
        mask = (sacc_bi == 0) & (blinks_bi == 0)
        #x_c = x_coord_sub[mask_test]
        x_coord_sub_m = x_coord_sub[mask]
        y_coord_sub_m = y_coord_sub[mask]
        #plt.plot(x_c, label='mask only sacc')
        #plt.plot(x_coord_sub_t, label='mask sacc + blinks')
        #plt.legend()
        #plt.title("Horizontal viewing behavior for two different masks")
        #plt.show()

        #print(x_coord_sub_m.shape)
        #print(y_coord_sub_m.shape)

        #plt.scatter(x_coord_sub, y_coord_sub, alpha=0.8, label='viewing')
        #plt.scatter(x_coord_sub_m, y_coord_sub_m, alpha= 0.5,label='masked viewing')
        #plt.title("Viewing behavior per subject per image")
        #plt.legend()
        #plt.show()


        print("Picture",sub_data['pic_#'])
        print("Category",sub_data['Category'])

        x_edg = np.linspace(192, 832, num=455)
        y_edg = np.linspace(144, 624, num=342)
        h_all_bins, xedg, yedg = np.histogram2d(x_coord_sub_m, y_coord_sub_m,  bins=[x_edg, y_edg])
        h_all = h_all_bins.T
        #H_flipped = np.flipud(H)
        H_smoothed = gaussian_filter(h_all, sigma=10)

        #print("X coordinate per subject per image", x_coord_sub)
        hist, xedges, yedges = np.histogram2d(x_coord_sub_m, y_coord_sub_m,  bins=[x_edges, y_edges])
        print("Hist", type(hist), hist.shape)
        print(hist)
        #print(hist.unique())
        H = hist.T
        print("Histogram flipped:", H)

        #get image that was shown for the averaged viewing behavior
        pic = sub_data['pic_#']
        pic = int(pic)
        c = sub_data['Category']
        for key, value in cat.items():
            if c == value:
                # if category of pic is matched in dictonary assign key (cat name) to first value
                c=key
        if c in ('landscapes', 'streets2'):
            img_path = os.path.join(pic_path, c, 'image'+str(pic)+'.png')
        else:
            img_path = os.path.join(pic_path, c, str(pic)+'.bmp')
        print(img_path)

        sub_list = [x for x in deep_gaze if c in x]
        #print(str(pic))
        deep_gaze_current_file = [x for x in sub_list if "_"+str(pic)+"." in x]
        print(deep_gaze_current_file[0], img_path)
        deep_gaze_path = os.path.join(path, deep_gaze_current_file[0])

        img = mpimg.imread(img_path)
        #deep_gaze_II = mpimg.imread(deep_gaze_path)
        deep_gaze_II = np.load(deep_gaze_path) #allow_pickle=True)
        print("DeepGaze", deep_gaze_II, deep_gaze_II.shape)

        #print("Check Nans", np.isnan(H_smoothed).any(), np.isnan(deep_gaze_II).any()) #--> no Nans
        #print("Check if STD is None or 0", H_smoothed.std(), deep_gaze_II.std())

        #deep_gaze_II = deep_gaze_II.T
        #reshape function doesn't work
        # numpy resize most likely not appropriate, check with Niels
        res = cv2.resize(deep_gaze_II, dsize=(3, 3), interpolation=cv2.INTER_AREA) # interpolation=cv2.INTER_AREA recommended for downsampling cv2.INTER_CUBIC
        print("Rescaled using opencv resize, based on interpolation", res.shape, res)
        res_deep = resize(deep_gaze_II, (3, 3))
        print("Rescaled using scikit-image resize", res_deep.shape, res_deep) #this might be the right one!
        #for displaying images I need to use nan, for correlation delete!
        res_deep[1,1] = np.nan
        #res[1,1] = np.nan
        #print(res_deep)
        #H_smoothed[1,1] = np.nan
        H[1,1] = np.nan
        #print("Histogram not smoothed", H)
        deep_8 = np.delete(res_deep, 4, None)
        print("Deep Gaze", deep_8, len(deep_8), type(deep_8))
        fdm = np.delete(H, 4, None)
        print("FDM", fdm, len(fdm), type(fdm))


        #fdm = H_smoothed.flatten()
        #print("Flattened FDM", fdm, type(fdm))
        corr_val, pearson_pval = stats.pearsonr(fdm, deep_8)#[0, 1] #[] gives me the value of the 2x2 array, so 0,1 is 0 rows and 1 column np.corrcoef
        # spearman correlation
        rho, pval = stats.spearmanr(fdm, deep_8)
        print("Correlation Pearson", corr_val)
        print("P-Value Pearson", pearson_pval)
        print("Correlation Spearman", rho)
        print("P-Value Spearman", pval)


        fig, ax = plt.subplots(2,2)
        ax[0][0].imshow(img, interpolation='nearest', extent=[xedg[0], xedg[-1],yedg[0],yedg[-1]]) #interpolation='gaussian'
        a0 = ax[0][0].imshow(H_smoothed, interpolation='nearest', origin = 'lower', extent=[xedg[0], xedg[-1],yedg[0],yedg[-1]], alpha=0.7) #cmap=plt.cm.RdBu
        ax[0][0].axis('off')
        divider1 = make_axes_locatable(ax[0][0])
        cax1 = divider1.append_axes("right", size="3%", pad=0.05)
        cbar1 = fig.colorbar(a0, ax=ax[0][0], cax = cax1)

        ax[0][1].imshow(img, interpolation='nearest', extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]]) #interpolation='gaussian'
        a1 = ax[0][1].imshow(H, interpolation='nearest', origin = 'lower', extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]], alpha=0.7) #cmap=plt.cm.RdBu
        ax[0][1].axis('off')
        divider2 = make_axes_locatable(ax[0][1])
        cax2 = divider2.append_axes("right", size="3%", pad=0.05)
        cbar2 = fig.colorbar(a1, ax=ax[0][1], cax = cax2)

        ax[1][0].imshow(img, interpolation='nearest', extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]]) #interpolation='gaussian'
        a2 = ax[1][0].imshow(deep_gaze_II, interpolation='nearest', origin = 'lower', extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]], alpha=0.7) #cmap=plt.cm.RdBu
        ax[1][0].axis('off')
        divider3 = make_axes_locatable(ax[1][0])
        cax3 = divider3.append_axes("right", size="3%", pad=0.05)
        cbar3 = fig.colorbar(a2, ax=ax[1][0], cax = cax3)

        ax[1][1].imshow(img, interpolation='nearest', extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]]) #interpolation='gaussian'
        a3 = ax[1][1].imshow(res_deep, interpolation='nearest', extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]], alpha=0.7) #cmap=plt.cm.RdBu
        ax[1][1].axis('off')
        divider4 = make_axes_locatable(ax[1][1])
        cax4 = divider4.append_axes("right", size="3%", pad=0.05)
        cbar4 = fig.colorbar(a3, ax=ax[1][1], cax = cax4)
        plt.show()






        dataset = pd.DataFrame({'FDM':fdm, 'Deep_gaze':deep_8, 'pic': sub_data['pic_#'], 'cat': sub_data['Category'], 'sub':name, 'pearson':corr_val, 'rho':rho}, columns=['FDM', 'Deep_gaze', "pic", "cat", "sub", "pearson", "rho"])
        print(dataset)

        #print("Zeros in deepgaze?", (dataset["Deep_gaze"] == 0).any()) #false --> no zeros in deep_gaze column Zeros in deepgaze? dataset["Deep_gaze"].isin([0]).any()
        #m1 = (dataset["FDM"] != 0).all() #false
        #m2 = (dataset["Deep_gaze"] != 0).all() #true
        #print(dataset["FDM"].iloc[dataset["FDM"].to_numpy().nonzero()[0]])
        #dataset_no_null = dataset[~(dataset == 0).any(axis=1)]
        #dataset_no_null = dataset[(dataset[['FDM','Deep_gaze']] != 0).all(axis=1)]
        #print(dataset_no_null.describe(), dataset_no_null_2.describe())
        #print(dataset_no_null.describe())
        #appended_data.append(dataset)
        #corr_values.append(corr_val)
        #rho_values.append(rho)



    #gives me dataframe over all images for one subject
    #df_final = pd.concat(appended_data) #ignore_index = True
    #print(df_final)
    #print(len(df_final['pic'].unique()),df_final['cat'].unique())
    #gives me list of corr values for all images within one subject
    #print("Corr_val mean within subject {}".format(name), np.mean(corr_values))
    #print("Rho value mean within subject {}".format(name), np.mean(rho_values))

    #pearson_mean.append(np.mean(corr_values))
    #rho_mean.append(np.mean(rho_values))
    #dataframes_sub.append(df_final)
#print(len(dataframes_sub))
#print(len(pearson_mean), len(rho_mean))
#print(pearson_mean)
#print(rho_mean)

#with open("rho_mean_vals.txt", "wb") as fp:
#    pickle.dump(rho_mean, fp)
#with open("pearson_mean_vals.txt", "wb") as fp:
#    pickle.dump(pearson_mean, fp)

#all_sub_all_img = pd.concat(dataframes_sub)
#print(all_sub_all_img) # save this
#all_sub_all_img.to_pickle('FDMs_all_sub')

        #dataset.hist(alpha=0.5, color='red')
        #dataset_no_null.hist(alpha=0.5,color='blue')
        #bins = np.linspace(0, 0.1, 100)
        #fig, ax = plt.subplots(1, 2)
        #ax = ax.ravel()

        #for idx in range(2):
        #    ax[idx].hist(dataset.iloc[:,idx], alpha=0.5, bins=100, density=True)
        #    ax[idx].hist(dataset_no_null.iloc[:,idx], alpha=0.3, bins=100, density=True) # bins=12
        #    ax[idx].set_title(dataset.columns[idx]+' with '+dataset_no_null.columns[idx]+" 0s removed")
            #ax[idx].legend(loc='upper left')
        #bins = np.linspace(0, 0.2, 100)

        #sns.distplot(dataset['FDM'], color="skyblue", label="With Zeros")
        #sns.distplot(dataset_no_null['FDM'], color="red", label="Zeros removed")
        #plt.legend()
        #plt.show()



        #if name in [53, 34, 24, 68, 35]:
        #plt.figure()
            #sns.lmplot(data=dataset, x='FDM', y='Deep_gaze')
            #plt.title(name)
            #plt.show()
            #plt.close()
        #corr_val = np.corrcoef(H_smoothed.flatten(), deep_gaze_II.flatten())[0, 1] #[] gives me the value of the 2x2 array, so 0,1 is 0 rows and 1 column
        #correlation_mat = ma.corrcoef(ma.masked_invalid(H_smoothed), ma.masked_invalid(deep_gaze_II))
        # spearman correlation
        #rho, pval = stats.spearmanr(H_smoothed.flatten(), deep_gaze_II.flatten())
        #print("Correlation Pearson", corr_val, corr_val.shape)
        #print("Correlation Spearman", rho, rho.shape)
        #print("P-Value Spearman", pval)
        #corr_values.append({
        #    'Subject': name,
        #    'Category': c,
        #    '#_pic':  pic,
        #    'pearson': corr_val,
        #    'rho': rho,
        #    'Image':sub_data['Images'] })

        #corr_values[(c,pic)] = corr_val

        #fig = plt.figure()
        #ax = fig.add_subplot(111)
        #cax = ax.matshow(corr_val, interpolation='nearest')
        #fig.colorbar(cax)
        #fig.suptitle('Correlation Matrix for img {} {}'.format(c, pic))
        #sns.heatmap(corr_val, xticklabels=corr_val.columns, yticklabels=corr_val.columns)#annot = True
        #plt.show()
        #plt.close()
        #corr_val = generate_correlation_map(H_smoothed, deep_gaze_II)
        #print("Correlation value of img {} {} for sub {}".format(c,pic, name), corr_val )
"""
"""
fig, axes = plt.subplots(nrows=2, ncols=2)
axes[0][0].imshow(img, interpolation='nearest', extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]]) #interpolation='gaussian'
axes[0][0].imshow(H_smoothed, interpolation='nearest', origin = 'lower', extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]], alpha=0.7, vmin=0, vmax=0.00001) #cmap=plt.cm.RdBu
axes[0][0].set_title('Viewing behavior of sub {} overlapped with image {} {}'.format(name, c, pic))
axes[0][0].axis('off')
#axes[1].set_xlabel('# of bins x axis')
#axes[1].set_ylabel('# of bins y axis')
# scaling of colorbar - 0-1.5 the best?
a0 = axes[0][1].imshow(H_smoothed, interpolation='nearest', origin = 'lower', extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]], vmin=0, vmax=0.00001)
divider1 = make_axes_locatable(axes[0][1])
cax1 = divider1.append_axes("right", size="3%", pad=0.05)
cbar1 = fig.colorbar(a0, ax=axes[0][1], cax = cax1)
axes[0][1].set_title('FDM of sub {} for image {} {}'.format(name, c, pic))
axes[0][1].axis('off')

#axes[0][2].imshow(img, interpolation='nearest') #interpolation='gaussian'
#axes[0][2].set_title('Original image {} {}'.format(pic_cat, name[1]))
#axes[0][2].axis('off')
#axes[2].set_title('imshow: interpolation gaussian for image # {}'.format(img))
#cbar1.ax.set_ylabel('Test', rotation = 270)

axes[1][0].imshow(img, interpolation='nearest', extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]]) #interpolation='gaussian'
axes[1][0].imshow(deep_gaze_II, interpolation='nearest', extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]], alpha=0.7) # cmap=plt.cm.RdBu  origin = 'lower',
axes[1][0].set_title('DeepGaze overlapped with image {} {}'.format(c, pic))
axes[1][0].axis('off')

a1 = axes[1][1].imshow(deep_gaze_II, interpolation='nearest', extent=[xedges[0], xedges[-1],yedges[0],yedges[-1]]) #origin = 'lower'
divider2 = make_axes_locatable(axes[1][1])
cax2 = divider2.append_axes("right", size="3%", pad=0.05)
cbar2 = fig.colorbar(a1, ax=axes[1][1], cax = cax2)
axes[1][1].set_title('DeepGaze density prediction for image {} {}'.format(c, pic))
axes[1][1].axis('off')

#axes[1][2].imshow(img, interpolation='nearest') #interpolation='gaussian'
#axes[1][2].set_title('Original image {} {}'.format(pic_cat, name[1]))
#axes[1][2].axis('off')
#cbar2.ax.set_ylabel('Test', rotation = 270)

plt.tight_layout()
#fig.colorbar(test)
#fig.colorbar(test, ax=axes.ravel().tolist(), shrink=0.5)

plt.show()
"""


#print("total number of runs, e.g. images or subjects", n) #150, so we have 150 groups, type df, of length 88 (subjects)
#test= data_image.apply(lambda x: x['Subjects'].unique())
# turn dictonary into df
#corr_val_df = pd.DataFrame(corr_values)
#print(corr_val_df)
#corr_val_df.to_pickle("corr_val_all.pkl")
#test = corr_val_df.drop(['Category', '#_pic', 'pval'], 1)
#test_group = test.groupby('Subject').mean().reset_index()
#print("Correlation mean per Subject", test_group)

#corr_val_sub = corr_val_df.groupby('Subject')
#test_1 = corr_val_df.drop(['Subject', 'pval'], 1)
#corr_val_img = corr_val_df.groupby(['Category', '#_pic']).mean().reset_index()
#print("Correlation mean per image", corr_val_img)
#for i, sub_data_df in corr_val_img:

#    pearson_mean = sub_data_df["pearson"].mean()
#    rho_mean = sub_data_df["rho"].mean()
#    print("Mean Correlation per Image {}".format(i), pearson_mean, rho_mean )
#for i, sub_data_df in corr_val_sub:
#    print("Subject", i)
#    pearson_mean = sub_data_df["pearson"].mean()
#    rho_mean = sub_data_df["rho"].mean()
#    print("Mean Correlation per Subject", pearson_mean, rho_mean )

#corr_val_df.hist(column='rho')
#test_group['pearson'].hist(legend=True)
#test_group['rho'].hist(alpha=0.5, legend=True)

#plt.title("Correlation averaged within subject")
#plt.legend()
#plt.show()


#corr_val_img['pearson'].hist(legend=True)
#corr_val_img['rho'].hist(alpha=0.5, legend=True)
#plt.xlabel('Subjects')
#plt.ylabel('Correlation')
#plt.title('Correlation averaged for image')
#plt.legend()
#plt.show()
#test = test.apply(pd.Series)
#print(type(test), test.head())
#corr_df = test.reset_index(inplace=True)
#print("Datafram Corr", type(corr_df), corr_df)

##test.apply(print)
#axarr = test_group.hist(column=['pearson', 'rho'], alpha=0.5, legend=True)

#for ax in axarr.flatten():
#    ax.set_xlabel("Subjects")
#    ax.set_ylabel("Correlation")
