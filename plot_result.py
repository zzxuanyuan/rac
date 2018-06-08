import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

''' block-size : b#, dict-chunk-size : c#, maxdict : d#, block-per-stripe : s#, compression level : e#, file name '''
COLUMN_NAMES=['block-size', 'dict-chunk-size', 'maxdict', 'block-per-stripe', 'compression-level', 'filename', 'filesize', 'dict-time', 'c-time-w-dict', 'seq-dec-time-w-dict', 'rand-dec-time-w-dict', 'c-time-wo-dict', 'seq-dec-time-wo-dict', 'rand-dec-time-wo-dict', 'rand-dec-size', 'c-size-w-dict','c-size-wo-dict']
df = pd.DataFrame(columns=COLUMN_NAMES)

''' level 1 with variant blocks per stripe '''
df.loc[0] = [4096, 4096, 4096, 8, 1, 'big', 6489102, 11.716381, 0.215606, 0.024390, 0.023741, 0.151077, 0.022626, 0.191955, 6488064, 4747232, 4120659]
df.loc[1] = [4096, 4096, 4096, 16, 1, 'big', 6489102, 12.596589, 0.239066, 0.027502, 0.027538, 0.143825, 0.024035, 0.384201, 6488064, 4619158, 3983998]
df.loc[2] = [4096, 4096, 4096, 32, 1, 'big', 6489102, 13.394346, 0.252131, 0.026614, 0.022910, 0.157389, 0.022340, 0.511773, 6488064, 4541398, 4001166]
df.loc[3] = [4096, 4096, 4096, 64, 1, 'big', 6489102, 13.762223, 0.239580, 0.029211, 0.027916, 0.157757, 0.018557, 1.126448, 6488064, 4496827, 3968896]
df.loc[4] = [4096, 4096, 4096, 128, 1, 'big', 6489102, 13.859363, 0.258491, 0.025811, 0.029535, 0.152448, 0.017326, 2.006949, 6488064, 4481268, 3954054]
df.loc[5] = [4096, 4096, 4096, 256, 1, 'big', 6489102, 14.975510, 0.240135, 0.024095, 0.024211, 0.150574, 0.016783, 4.080633, 6488064, 4486367, 3946584]

''' level 1 with variant block size '''
df.loc[6] = [8192, 4096, 4096, 8, 1, 'big', 6489102, 12.569408, 0.213014, 0.029521, 0.024685, 0.135421, 0.023684, 0.180112, 6488064, 4509434, 3983998]
df.loc[7] = [16384, 4096, 4096, 8, 1, 'big', 6489102, 12.388515, 0.186539, 0.021641, 0.018246, 0.166875, 0.016896, 0.138512, 6488064, 4331026, 4001166]
df.loc[8] = [32768, 4096, 4096, 8, 1, 'big', 6489102, 13.398630, 0.178143, 0.024669, 0.019897, 0.172245, 0.016391, 0.129033, 6488064, 4172934, 3968896]
df.loc[9] = [65536, 4096, 4096, 8, 1, 'big', 6489102, 14.180327, 0.179243, 0.020662, 0.022820, 0.159134, 0.017278, 0.138316, 6488064, 4065039, 3954054]
df.loc[10] = [131072, 4096, 4096, 8, 1, 'big', 6489102, 14.936123, 0.160111, 0.021032, 0.017550, 0.151350, 0.016510, 0.123637, 6422528, 4005005, 3946584]

''' level 1 with variant dictionary size '''
df.loc[11] = [4096, 4096, 8192, 256, 1, 'big', 6489102, 14.454033, 0.277840, 0.022033, 0.019068, 0.167079, 0.016763, 4.055762, 6488064, 4303212, 3946584]
df.loc[12] = [4096, 4096, 16384, 256, 1, 'big', 6489102, 14.144350, 0.372117, 0.026970, 0.023551, 0.154270, 0.017535, 3.972519, 6488064, 4172853, 3946584]
df.loc[13] = [4096, 4096, 32768, 256, 1, 'big', 6489102, 14.063793, 0.558558, 0.023274, 0.019313, 0.151908, 0.017130, 4.068578, 6488064, 4164880, 3946584]
df.loc[14] = [4096, 4096, 65536, 256, 1, 'big', 6489102, 13.651855, 0.910053, 0.028467, 0.018744, 0.156575, 0.016495, 4.051256, 6488064, 4228438, 3946584]
df.loc[15] = [4096, 4096, 131072, 256, 1, 'big', 6489102, 12.988945, 0.880230, 0.033055, 0.033543, 0.168966, 0.019042, 4.032770, 6488064, 4711450, 3946584]

''' level 1 with variant dictionary chunk size '''
df.loc[16] = [4096, 8192, 65536, 256, 1, 'big', 6489102, 13.666903, 0.881506, 0.028884, 0.025097, 0.173440, 0.018862, 3.994764, 6488064, 4228708, 3946584]
df.loc[17] = [4096, 16384, 65536, 256, 1, 'big', 6489102, 13.443700, 0.870179, 0.026648, 0.017719, 0.151370, 0.016701, 4.052547, 6488064, 4230164, 3946584]
df.loc[18] = [4096, 32768, 65536, 256, 1, 'big', 6489102, 12.977094, 0.874804, 0.032881, 0.023995, 0.159671, 0.018625, 4.002807, 6488064, 4227039, 3946584]

''' level 3 with variant blocks per stripe '''
df.loc[19] = [4096, 4096, 4096, 8, 3, 'big', 6489102, 13.541882, 0.631895, 0.039732, 0.033864, 0.490164, 0.022073, 0.170340, 6488064, 4126363, 3505686]
df.loc[20] = [4096, 4096, 4096, 16, 3, 'big', 6489102, 13.548690, 0.633929, 0.039195, 0.018161, 0.476446, 0.022870, 0.316772, 6488064, 3949805, 3247727]
df.loc[21] = [4096, 4096, 4096, 32, 3, 'big', 6489102, 14.109613, 0.652278, 0.040661, 0.038934, 0.478794, 0.019828, 0.605680, 6488064, 3847837, 3085921]
df.loc[22] = [4096, 4096, 4096, 64, 3, 'big', 6489102, 14.679614, 0.651037, 0.028836, 0.037750, 0.485236, 0.022208, 1.273272, 6488064, 3801214, 3000443]
df.loc[23] = [4096, 4096, 4096, 128, 3, 'big', 6489102, 15.369518, 0.683444, 0.028799, 0.035938, 0.517798, 0.021701, 2.388095, 6488064, 3783550, 2961296]
df.loc[24] = [4096, 4096, 4096, 256, 3, 'big', 6489102, 16.528321, 0.663147, 0.036483, 0.041030, 0.524071, 0.029541, 4.867560, 6488064, 3792318, 2940710]

''' level 3 with variant block size '''
df.loc[25] = [8192, 4096, 4096, 8, 3, 'big', 6489102, 14.103217, 0.649336, 0.033794, 0.037322, 0.456928, 0.020123, 0.146991, 6488064, 3809155, 3247727]
df.loc[26] = [16384, 4096, 4096, 8, 3, 'big', 6489102, 13.908722, 0.692411, 0.031008, 0.030626, 0.473853, 0.025337, 0.165353, 6488064, 3532048, 3085921]
df.loc[27] = [32768, 4096, 4096, 8, 3, 'big', 6489102, 15.090300, 0.896262, 0.024569, 0.030162, 0.488841, 0.021410, 0.155526, 6488064, 3274910, 3000443]
df.loc[28] = [65536, 4096, 4096, 8, 3, 'big', 6489102, 15.107022, 1.179296, 0.027373, 0.026245, 0.501949, 0.025190, 0.168207, 6488064, 3037309, 2961296]
df.loc[29] = [131072, 4096, 4096, 8, 3, 'big', 6489102, 16.285958, 1.306184, 0.027039, 0.027789, 0.475875, 0.019234, 0.159924, 6488064, 2889204, 2940710]

''' level 3 with variant dictionary size '''
df.loc[30] = [4096, 4096, 8192, 256, 3, 'big', 6489102, 15.466024, 0.908078, 0.027682, 0.028966, 0.478925, 0.023261, 4.728335, 6488064, 3531702, 2940710]
df.loc[31] = [4096, 4096, 16384, 256, 3, 'big', 6489102, 15.651928, 1.348537, 0.032140, 0.022869, 0.483059, 0.019297, 4.996209, 6488064, 3291333, 2940710]
df.loc[32] = [4096, 4096, 32768, 256, 3, 'big', 6489102, 15.097456, 2.024780, 0.031656, 0.024088, 0.481447, 0.019354, 4.802181, 6488064, 3102127, 2940710]
df.loc[33] = [4096, 4096, 65536, 256, 3, 'big', 6489102, 14.824586, 3.079823, 0.025696, 0.021882, 0.489474, 0.019370, 4.760722, 6488064, 3034525, 2940710]
df.loc[34] = [4096, 4096, 131072, 256, 3, 'big', 6489102, 14.161229, 3.066829, 0.031941, 0.039027, 0.484678, 0.020881, 4.761188, 6488064, 3544297, 2940710]
#df.loc[5] = [4096, 4096, 262144, 256, 3, 'big', 6489102, 14.196185, 3.071249, 0.038304, 0.037363, 0.484556, 0.025371, 4.898985, 6488064, 4527397, 2940710]

''' level 3 with variant dictionary chunk size '''
df.loc[35] = [4096, 8192, 65536, 256, 3, 'big', 6489102, 14.612903, 3.139835, 0.029300, 0.021663, 0.484277, 0.019673, 4.784113, 6488064, 3034638, 2940710]
df.loc[36] = [4096, 16384, 65536, 256, 3, 'big', 6489102, 14.215566, 3.091051, 0.029621, 0.019961, 0.495184, 0.019847, 4.770087, 6488064, 3035729, 2940710]
df.loc[37] = [4096, 32768, 65536, 256, 3, 'big', 6489102, 14.203804, 3.133730, 0.028117, 0.028923, 0.481301, 0.023324, 4.733170, 6488064, 3036873, 2940710]

''' osdb ''' 
''' level 1 with variant blocks per stripe '''
df.loc[38] = [4096, 4096, 4096, 8, 1, 'osdb', 10085684, 18.701072, 0.328927, 0.027485, 0.026144, 0.200143, 0.024127, 0.185470, 10084352, 7375072, 5884245]
df.loc[39] = [4096, 4096, 4096, 16, 1, 'osdb', 10085684, 19.137210, 0.341186, 0.030018, 0.025080, 0.188874, 0.029001, 0.341340, 10084352, 7003371, 5246156]
df.loc[40] = [4096, 4096, 4096, 32, 1, 'osdb', 10085684, 19.938976, 0.347508, 0.027959, 0.032532, 0.224875, 0.026218, 0.629424, 10084352, 6803023, 5438280]
df.loc[41] = [4096, 4096, 4096, 64, 1, 'osdb', 10085684, 20.811047, 0.349873, 0.026016, 0.025010, 0.216521, 0.029646, 1.226797, 10084352, 6694676, 5345984]
df.loc[42] = [4096, 4096, 4096, 128, 1, 'osdb', 10085684, 22.411599, 0.345139, 0.030088, 0.027273, 0.220418, 0.022252, 2.538609, 10084352, 6627898, 5301712]
df.loc[43] = [4096, 4096, 4096, 256, 1, 'osdb', 10085684, 23.931411, 0.360889, 0.030510, 0.027859, 0.209585, 0.026527, 4.893201, 10084352, 6591543, 5278054]

''' level 1 with variant block size '''
df.loc[44] = [8192, 4096, 4096, 8, 1, 'osdb', 10085684, 19.866215, 0.315463, 0.034790, 0.024402, 0.186344, 0.024052, 0.171219, 10084352, 6691133, 5246156]
df.loc[45] = [16384, 4096, 4096, 8, 1, 'osdb', 10085684, 19.894911, 0.274029, 0.030589, 0.026082, 0.211239, 0.026842, 0.158851, 10076160, 6078830, 5438280]
df.loc[46] = [32768, 4096, 4096, 8, 1, 'osdb', 10085684, 20.803073, 0.240677, 0.031450, 0.020603, 0.214844, 0.021429, 0.175071, 10059776, 5654060, 5345984]
df.loc[47] = [65536, 4096, 4096, 8, 1, 'osdb', 10085684, 22.245988, 0.225438, 0.025599, 0.020647, 0.201641, 0.024135, 0.157699, 10027008, 5445327, 5301712]
df.loc[48] = [131072, 4096, 4096, 8, 1, 'osdb', 10085684, 23.531033, 0.225700, 0.027067, 0.028891, 0.213290, 0.021042, 0.173125, 9961472, 5346990, 5278054]

''' level 1 with variant dictionary size '''
df.loc[49] = [4096, 4096, 8192, 256, 1, 'osdb', 10085684, 22.985590, 0.409239, 0.031261, 0.030170, 0.214488, 0.021048, 4.883159, 10084352, 5591691, 5278054]
df.loc[50] = [4096, 4096, 16384, 256, 1, 'osdb', 10085684, 22.284024, 0.507390, 0.033380, 0.023826, 0.214416, 0.022735, 5.214614, 10084352, 4969980, 5278054]
df.loc[51] = [4096, 4096, 32768, 256, 1, 'osdb', 10085684, 21.263986, 0.784359, 0.029996, 0.027076, 0.216242, 0.026416, 5.023053, 10084352, 4882574, 5278054]
df.loc[52] = [4096, 4096, 65536, 256, 1, 'osdb', 10085684, 21.455761, 1.383469, 0.029839, 0.025946, 0.210055, 0.027925, 4.836884, 10084352, 5044500, 5278054]
df.loc[53] = [4096, 4096, 131072, 256, 1, 'osdb', 10085684, 20.605187, 1.299979, 0.029996, 0.025987, 0.224349, 0.027132, 4.884519, 10084352, 5719136, 5278054]

''' level 1 with variant dictionary chunk size '''
df.loc[54] = [4096, 8192, 65536, 256, 1, 'osdb', 10085684, 20.846204, 1.335272, 0.028739, 0.026224, 0.225984, 0.021833, 5.016333, 10084352, 5071080, 5278054]
df.loc[55] = [4096, 16384, 65536, 256, 1, 'osdb', 10085684, 20.481975, 1.320190, 0.030447, 0.021987, 0.213685, 0.026524, 4.975645, 10084352, 5041907, 5278054]
df.loc[56] = [4096, 32768, 65536, 256, 1, 'osdb', 10085684, 20.211560, 1.291746, 0.031653, 0.024320, 0.219180, 0.022440, 4.868705, 10084352, 5050899, 5278054]

''' level 3 with variant blocks per stripe '''
'''
df.loc[19] = [4096, 4096, 4096, 8, 3, 10085684]
df.loc[20] = [4096, 4096, 4096, 16, 3, 10085684]
df.loc[21] = [4096, 4096, 4096, 32, 3, 10085684]
df.loc[22] = [4096, 4096, 4096, 64, 3, 10085684]
df.loc[23] = [4096, 4096, 4096, 128, 3, 10085684]
df.loc[24] = [4096, 4096, 4096, 256, 3, 10085684]
'''
''' level 3 with variant block size '''
'''
df.loc[25] = [8192, 4096, 4096, 8, 3, 10085684]
df.loc[26] = [16384, 4096, 4096, 8, 3, 10085684]
df.loc[27] = [32768, 4096, 4096, 8, 3, 10085684]
df.loc[28] = [65536, 4096, 4096, 8, 3, 10085684]
df.loc[29] = [131072, 4096, 4096, 8, 3, 10085684]
'''
''' level 3 with variant dictionary size '''
'''
df.loc[30] = [4096, 4096, 8192, 256, 3, 10085684]
df.loc[31] = [4096, 4096, 16384, 256, 3, 10085684]
df.loc[32] = [4096, 4096, 32768, 256, 3, 10085684]
df.loc[33] = [4096, 4096, 65536, 256, 3, 10085684]
df.loc[34] = [4096, 4096, 131072, 256, 3, 10085684]
'''
''' level 3 with variant dictionary chunk size '''
'''
df.loc[35] = [4096, 8192, 65536, 256, 3, 10085684]
df.loc[36] = [4096, 16384, 65536, 256, 3, 10085684]
df.loc[37] = [4096, 32768, 65536, 256, 3, 10085684]
'''
''' ooffice ''' 
''' level 1 with variant blocks per stripe '''
df.loc[57] = [4096, 4096, 4096, 8, 1, 'ooffice', 6152192, 11.580942, 0.193354, 0.017001, 0.014413, 0.165375, 0.020028, 0.150876, 6152192, 4757679, 4263894]
df.loc[58] = [4096, 4096, 4096, 16, 1, 'ooffice', 6152192, 11.809262, 0.207248, 0.016542, 0.016464, 0.151168, 0.017016, 0.258256, 6152192, 4665060, 4202445]
df.loc[59] = [4096, 4096, 4096, 32, 1, 'ooffice', 6152192, 12.245707, 0.208888, 0.018604, 0.014717, 0.164823, 0.013002, 0.402193, 6152192, 4603416, 4364392]
df.loc[60] = [4096, 4096, 4096, 64, 1, 'ooffice', 6152192, 12.317984, 0.207309, 0.020427, 0.014744, 0.163688, 0.014952, 0.803677, 6152192, 4588625, 4354814]
df.loc[61] = [4096, 4096, 4096, 128, 1, 'ooffice', 6152192, 13.102338, 0.229263, 0.017112, 0.014946, 0.163386, 0.013472, 1.581143, 6152192, 4583796, 4345650]
df.loc[62] = [4096, 4096, 4096, 256, 1, 'ooffice', 6152192, 13.835337, 0.226176, 0.019723, 0.014914, 0.168355, 0.013154, 3.412196, 6152192, 4585952, 4341223]

''' level 1 with variant block size '''
df.loc[63] = [8192, 4096, 4096, 8, 1, 'ooffice', 6152192, 11.630635, 0.184409, 0.020001, 0.015738, 0.148318, 0.016253, 0.137248, 6152192, 4604906, 4202445]
df.loc[64] = [16384, 4096, 4096, 8, 1, 'ooffice', 6152192, 12.081106, 0.179548, 0.016111, 0.014318, 0.187261, 0.013745, 0.120115, 6144000, 4498057, 4364392]
df.loc[65] = [32768, 4096, 4096, 8, 1, 'ooffice', 6152192, 13.036529, 0.174082, 0.014951, 0.014066, 0.153957, 0.012981, 0.101018, 6127616, 4445161, 4354814]
df.loc[66] = [65536, 4096, 4096, 8, 1, 'ooffice', 6152192, 13.051825, 0.156458, 0.020142, 0.014201, 0.172035, 0.017139, 0.100043, 6094848, 4406543, 4345650]
df.loc[67] = [131072, 4096, 4096, 8, 1, 'ooffice', 6152192, 13.887253, 0.160815, 0.016899, 0.015878, 0.153406, 0.014033, 0.107319, 6029312, 4376205, 4341223]

''' level 1 with variant dictionary size '''
df.loc[68] = [4096, 4096, 8192, 256, 1, 'ooffice', 6152192, 13.234216, 0.279050, 0.018035, 0.014172, 0.160628, 0.013094, 3.213210, 6152192, 4543989, 4341223]
df.loc[69] = [4096, 4096, 16384, 256, 1, 'ooffice', 6152192, 12.879771, 0.353325, 0.021015, 0.015579, 0.154642, 0.013405, 3.308654, 6152192, 4506971, 4341223]
df.loc[70] = [4096, 4096, 32768, 256, 1, 'ooffice', 6152192, 12.598406, 0.450798, 0.018931, 0.020134, 0.157144, 0.013763, 3.323311, 6152192, 4938737, 4341223]
df.loc[71] = [4096, 4096, 65536, 256, 1, 'ooffice', 6152192, 12.186997, 0.849932, 0.018724, 0.016465, 0.157762, 0.012994, 3.183434, 6152192, 4624643, 4341223]
df.loc[72] = [4096, 4096, 131072, 256, 1, 'ooffice', 6152192, 13.807419, 0.917171, 0.022583, 0.022351, 0.181987, 0.018031, 3.314355, 6152192, 5458972, 4341223]

''' level 1 with variant dictionary chunk size '''
df.loc[73] = [4096, 8192, 65536, 256, 1, 'ooffice', 6152192, 11.850398, 0.849755, 0.018027, 0.016142, 0.158558, 0.014495, 3.271279, 6152192, 4614155, 4341223]
df.loc[74] = [4096, 16384, 65536, 256, 1, 'ooffice', 6152192, 11.601792, 0.715383, 0.018545, 0.015156, 0.158426, 0.014421, 3.154924, 6152192, 5038686, 4341223]
df.loc[75] = [4096, 32768, 65536, 256, 1, 'ooffice', 6152192, 11.732358, 0.726204, 0.018614, 0.018474, 0.157823, 0.015045, 3.218587, 6152192, 5041325, 4341223]

''' level 3 with variant blocks per stripe '''
'''
df.loc[19] = [4096, 4096, 4096, 8, 3, ]
df.loc[20] = [4096, 4096, 4096, 16, 3, ]
df.loc[21] = [4096, 4096, 4096, 32, 3, ]
df.loc[22] = [4096, 4096, 4096, 64, 3, ]
df.loc[23] = [4096, 4096, 4096, 128, 3, ]
df.loc[24] = [4096, 4096, 4096, 256, 3, ]
'''
''' level 3 with variant block size '''
'''
df.loc[25] = [8192, 4096, 4096, 8, 3, ]
df.loc[26] = [16384, 4096, 4096, 8, 3, ]
df.loc[27] = [32768, 4096, 4096, 8, 3, ]
df.loc[28] = [65536, 4096, 4096, 8, 3, ]
df.loc[29] = [131072, 4096, 4096, 8, 3, ]
'''
''' level 3 with variant dictionary size '''
'''
df.loc[30] = [4096, 4096, 8192, 256, 3, ]
df.loc[31] = [4096, 4096, 16384, 256, 3, ]
df.loc[32] = [4096, 4096, 32768, 256, 3, ]
df.loc[33] = [4096, 4096, 65536, 256, 3, ]
df.loc[34] = [4096, 4096, 131072, 256, 3, ]
'''
''' level 3 with variant dictionary chunk size '''
'''
df.loc[35] = [4096, 8192, 65536, 256, 3, ]
df.loc[36] = [4096, 16384, 65536, 256, 3, ]
df.loc[37] = [4096, 32768, 65536, 256, 3, ]
'''

def autolabel(rects, ax):
	"""
	Attach a text label above each bar displaying its height
	"""
	for rect in rects:
		height = rect.get_height()
		ax.text(rect.get_x() + rect.get_width()/2., 1.05*height, '%d' % int(height), ha='center', va='bottom')

def plot_comp(seriesD, seriesX, seriesY, xLabel, yLabel, title, ticklabels, legend, log=0):
	N = seriesX.shape[0]
	print "size is ", N
	ind = np.arange(N)  # the x locations for the groups
	width = 0.35       # the width of the bars

	fig, ax = plt.subplots()
	rects1 = ax.bar(ind, seriesX, width, color='orange', label=legend[0])
	rects2 = ax.bar(ind, seriesD, width, bottom=seriesX, color='g', label=legend[1])
	rects3 = ax.bar(ind + width, seriesY, width, color='b', label=legend[2])

	# add some text for labels, title and axes ticks
	ax.set_xlabel(xLabel, fontsize=14)
	ax.set_ylabel(yLabel, fontsize=14)
#	ax.set_title(title)
	ax.set_xticks(ind + width/2)
	ax.set_xticklabels(ticklabels, fontsize=14)
	if log == 1:
		ax.set_yscale("log")
	plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0., fontsize=12)
	plt.show()


def plot_dec(seriesX, seriesY, xLabel, yLabel, title, ticklabels, legend, log=0):
	N = seriesX.shape[0]
	print "size is ", N
	ind = np.arange(N)  # the x locations for the groups
	width = 0.35       # the width of the bars

	fig, ax = plt.subplots()
	rects1 = ax.bar(ind, seriesX, width, color='orange', label=legend[0])
	rects2 = ax.bar(ind + width, seriesY, width, color='b', label=legend[1])

	# add some text for labels, title and axes ticks
	ax.set_xlabel(xLabel, fontsize=14)
	ax.set_ylabel(yLabel, fontsize=14)
	ax.set_xticks(ind + width/2)
	ax.set_xticklabels(ticklabels, fontsize=14)
	if log == 1:
		ax.set_yscale("log")
	plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0., fontsize=12)
	plt.show()

df_c1 = df[(df['block-size']==4096) & (df['dict-chunk-size']==4096) & (df['maxdict']==4096) & (df['compression-level']==1) & (df['filename']=='osdb')]
print df_c1['filename']
df_c3 = df[(df['block-size']==4096) & (df['dict-chunk-size']==4096) & (df['maxdict']==4096) & (df['compression-level']==3)]

plot_comp(df_c1['dict-time'], df_c1['c-time-w-dict'], df_c1['c-time-wo-dict'], 'Number of Blocks per Stripe', 'Time to Compress Data (Second)', 'Variant Blocks-Per-Stripe', ('8','16','32','64','128','256'), ('rac fast', 'dictionary fast', 'normal fast'))
plot_dec(df_c1['seq-dec-time-w-dict'], df_c1['seq-dec-time-wo-dict'], 'Number of Blocks per Stripe', 'Time to Sequentially Decompress Data (Second)', 'Variant Blocks-Per-Stripe', ('8','16','32','64','128','256'), ('rac fast', 'normal fast'))
plot_dec(df_c1['rand-dec-time-w-dict'], df_c1['rand-dec-time-wo-dict'], 'Number of Blocks per Stripe', 'Time to Randomly Decompress Data (Second)', 'Variant Blocks-Per-Stripe', ('8','16','32','64','128','256'), ('rac fast', 'normal fast'))
plot_dec(df_c1['filesize']/df_c1['c-size-w-dict'], df_c1['filesize']/df_c1['c-size-wo-dict'], 'Number of Blocks per Stripe', 'Compression Ratio', 'Variant Blocks-Per-Stripe', ('8','16','32','64','128','256'), ('rac fast', 'normal fast'))

#plot_dec(df_c3['filesize']/df_c3['c-size-w-dict'], df_c3['filesize']/df_c3['c-size-wo-dict'], 'Number of Blocks per Stripe', 'Compression Ratio', 'Variant Blocks-Per-Stripe', ('8','16','32','64','128','256'), ('rac hc', 'normal hc'))

df_c1 = df[(df['block-per-stripe']==8) & (df['dict-chunk-size']==4096) & (df['maxdict']==4096) & (df['compression-level']==1) & (df['filename']=='osdb')]
df_c3 = df[(df['block-per-stripe']==8) & (df['dict-chunk-size']==4096) & (df['maxdict']==4096) & (df['compression-level']==3)]

plot_comp(df_c1['dict-time'], df_c1['c-time-w-dict'], df_c1['c-time-wo-dict'],  'Block Size', 'Time to Compress Data (Second)', 'Variant Block-Size', ('4K','8K','16K','32K','64K','128K'), ('rac fast', 'dictionary fast', 'normal fast'))
plot_dec(df_c1['seq-dec-time-w-dict'], df_c1['seq-dec-time-wo-dict'], 'Block Size', 'Time to Sequentially Decompress Data (Second)', 'Variant Block-Size', ('4K','8K','16K','32K','64K','128K'), ('rac fast', 'normal fast'))
plot_dec(df_c1['rand-dec-time-w-dict'], df_c1['rand-dec-time-wo-dict'], 'Block Size', 'Time to Randomly Decompress Data (Second)', 'Variant Block-Size', ('4K','8K','16K','32K','64K','128K'), ('rac fast', 'normal fast'))
plot_dec(df_c1['filesize']/df_c1['c-size-w-dict'], df_c1['filesize']/df_c1['c-size-wo-dict'], 'Block Size', 'Compression Ratio', 'Variant Blocks-Size', ('4K', '8K','16K','32K','64K','128K'), ('rac fast', 'normal fast'))

#plot_dec(df_c3['filesize']/df_c3['c-size-w-dict'], df_c3['filesize']/df_c3['c-size-wo-dict'], 'Block Size', 'Compression Ratio', 'Variant Blocks-Per-Stripe', ('4K', '8K','16K','32K','64K','128K'), ('rac hc', 'normal hc'))
df_c1 = df[(df['block-per-stripe']==256) & (df['dict-chunk-size']==4096) & (df['block-size']==4096) & (df['compression-level']==1) & (df['filename']=='osdb')]
df_c3 = df[(df['block-per-stripe']==256) & (df['dict-chunk-size']==4096) & (df['block-size']==4096) & (df['compression-level']==3)]

plot_comp(df_c1['dict-time'], df_c1['c-time-w-dict'], df_c1['c-time-wo-dict'],  'Dictionary Size', 'Time to Compress Data (Second)', 'Variant Dictionary-Size', ('4K','8K','16K','32K','64K','128K'), ('rac fast', 'dictionary fast', 'normal fast'))
plot_dec(df_c1['seq-dec-time-w-dict'], df_c1['seq-dec-time-wo-dict'], 'Dictionary Size', 'Time to Sequentially Decompress Data (Second)', 'Variant Dictionary-Size', ('4K','8K','16K','32K','64K','128K'), ('rac fast', 'normal fast'))
plot_dec(df_c1['rand-dec-time-w-dict'], df_c1['rand-dec-time-wo-dict'], 'Dictionary Size', 'Time to Randomly Decompress Data (Second)', 'Variant Dictionary-Size', ('4K','8K','16K','32K','64K','128K'), ('rac fast', 'normal fast'))
plot_dec(df_c1['filesize']/df_c1['c-size-w-dict'], df_c1['filesize']/df_c1['c-size-wo-dict'], 'Dictionary Size', 'Compression Ratio', 'Variant Dictionary-Size', ('4K', '8K','16K','32K','64K','128K'), ('rac fast', 'normal fast'))