# -*-coding:UTF-8-*-
# Author  : xuzhou
# Date    : 2019-07-15
# Function: 计算砖和砂浆的CDP模型材料参数。
# 运行方式：1)CAE界面内，菜单栏File->Run Script...
#          2)启动Abaqus Command,输入: abaqus cae noGUI=E:/...../fileName.py，注意CMD工作目录为结果文件存储目录
#===================================================================================================
from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
import displayGroupOdbToolset as dgo
import customKernel
import xyPlot
# ==========================================【初始化】==============================================
import numpy as np
import os,sys,re,time
# ==================================================================================================
outputDir = "D:/temp/Wall-test-xuzhou/W1/outputInp"         #输出的inp文件目录
caeFullPath = "D:/temp/Wall-test-xuzhou/W1/W1-0616.cae"     #获得cae的完整路径
matName_brick = "brick"                                     #砖材料名称（所有cae文件中的材料名称都要一致）
matName_mortar = "mortar"                                   #水泥材料名称（所有cae文件中的材料名称都要一致）
#matName_uhpc = "UHPC"                                      #UHPC材料名称（所有cae文件中的材料名称都要一致）

#从csv文件中读取工况数据
fileName = "C:/Users/ASUS/Desktop/fc_Brick_Lime.csv"          #指定文件路径
fc_Brick_Lime = []                                            #初始化数组
f = open(fileName,'r')                                        #打开文件
fstr = f.readlines()                                          #读取文件存储到fstr
f.close()                                                     #文件关闭
for line in fstr[1:]:                                         #循环文件中各行，排除第一行
    words = line.replace("\n","").split(",")                  #擅长换行符，并以逗号分隔数据成列表
    if len(words)==7:                                         #当列表长度为9时，表示有数据
        fc_Brick_Lime.append([float(words[0]),float(words[1]),float(words[2]),float(words[3]),float(words[4]),
        float(words[5]),float(words[6])])                     #存储数据到列表

print fc_Brick_Lime             #测试打印数据

#砖材料参数输入
den_brick = 1800                #密度
elas_brick = 6.673e9            #弹性模量
pos_brick = 0.15                #泊松比
alpha_a_brick = 2.2             #砖压缩公式中的alpha_a
alpha_d_brick = 2               #砖压缩公式中的alpha_d
ec0_brick = 0.006               #极限压应力对应的应变
ecp_brick = 0.03                #取值最大压应变
et0_brick = 0.0001              #极限拉应力对应的应变
etp_brick = 0.001               #取值最大拉应变
p1_brick = 30                   #CDP模型中的 膨胀角
p2_brick = 0.1                  #CDP模型中的 偏心率
p3_brick = 1.16                 #CDP模型中的 双向/单向抗压强度比
p4_brick = 0.6667               #CDP模型中的 不变应力比
p5_brick = 0.005                #CDP模型中的 粘聚系数
wt_brick = 0.6                  #拉伸恢复系数
wc_brick = 0                    #压缩恢复系数
fcInitial_brick = 0.08          #初始损伤的压应力比率
ftInitial_brick = 0.95          #初始损伤的拉应力比率
dfct_brick = 0.05               #单个增量步运行的压应力变化比例（控制输出点间距）
maxDamage_brick = 0.9           #最大损伤阈值

#砂浆材料参数输入
den_mortar = 2000               #密度
elas_mortar = 100e6             #弹性模量
pos_mortar = 0.2                #泊松比
fc_mortar = 0.25e6              #抗压强度
ft_mortar = 0.025e6             #抗拉强度
ec0_mortar = 0.008              #极限压应力对应的应变
ecp_mortar = 0.05               #取值最大压应变
p1_mortar = 30                  #CDP模型中的 膨胀角
p2_mortar = 0.1                 #CDP模型中的 偏心率
p3_mortar = 1.16                #CDP模型中的 双向/单向抗压强度比
p4_mortar = 0.6667              #CDP模型中的 不变应力比
p5_mortar = 0.005               #CDP模型中的 粘聚系数
wt_mortar = 1                   #拉伸恢复系数
wc_mortar = 0                   #压缩恢复系数
fcInitial_mortar = 0.2          #初始损伤的压应力比率
dfct_mortar = 0.05              #单个增量步运行的压应力变化比例（控制输出点间距）
maxDamage_mortar = 0.9          #最大损伤阈值


class xuzhouCDP_Brick():
    def __init__( self ):
        self.xc=[]              #参数记录压应变
        self.xt=[]              #参数记录拉应变
        self.yc=[]              #参数记录压应力
        self.yt=[]              #参数记录拉应力
        self.dc=[]              #参数记录压损伤
        self.dt=[]              #参数记录拉损伤
        self.spc=[]             #参数记录压缩残余应变
        self.spt=[]             #参数记录拉伸残余应变
        
    def createMaterial(self,matName,den,elas,pos,fc,ft,p1,p2,p3,p4,p5,wt,wc,
        fcInitial,ftInitial,dfct,maxDamage,alpha_a,alpha_d,ec0,ecp,et0,etp,modelName):
        self.matName = matName          #材料名称
        self.den = den                  #密度
        self.elas = elas                #弹性模量
        self.pos = pos                  #泊松比
        self.fc = fc                    #压缩极限应力
        self.fcInitial = fcInitial      #初始压缩屈服应力
        self.ft = ft                    #拉伸极限应力
        self.ftInitial = ftInitial      #初始拉伸屈服应力
        self.alpha_a = alpha_a          #砖压缩公式中的alpha_a
        self.alpha_d = alpha_d          #砖压缩公式中的alpha_d
        self.ec0 = ec0                  #极限压应力对应的应变
        self.ecp = ecp                  #取值最大压应变
        self.et0 = et0                  #极限拉应力对应的应变
        self.etp = etp                  #取值最大拉应变
        self.p1 = p1                    #膨胀角
        self.p2 = p2                    #偏心系数
        self.p3 = p3                    #双向/单向抗压强度比
        self.p4 = p4                    #不变应力比
        self.p5 = p5                    #粘聚系数
        self.wt = wt                    #拉伸恢复系数
        self.wc = wc                    #压缩恢复系数
        self.dfct = dfct                #用于控制输出，当前后两个数据变化大于该数据时设置输出
        self.maxDamage = maxDamage      #允许的最大损伤
        self.modelName = modelName      #材料定义的模型名称
        
        #【计算塑性】
        if self.cal_Plastic()==0:
            return 0
        
        #【计算损伤】
        self.cal_DamageEnergy()
        
        #【数据截断】
        self.yc = [round(i,5) for i in self.yc]
        self.yt = [round(i,5) for i in self.yt]
        self.dc = [round(i,8) for i in self.dc]
        self.dt = [round(i,8) for i in self.dt]
        self.inpc = [round(i,8) for i in self.inpc]
        self.inpt = [round(i,8) for i in self.inpt]
        self.elas=self.elas
        #【生成材料】
        self.createAbaqusMat()
        print "Finish"
        return 1
    
    def cal_Plastic(self):
        #计算砖的塑性
        #========================压应力-应变曲线========================
        xc = [0]
        yc = [0]
        dx = 0                                                          #前期线性段引起应力应变曲线左移量
        x1 = np.array(range(1001))/1000.0                               #归一化应变
        y = self.alpha_a*x1+(3-2*self.alpha_a)*x1**2+(self.alpha_a-2)*x1**3 #归一化应力
        checkElas_c = self.ec0/(self.fc/self.elas)                      #初始模量的归一化变换
        for i in range(len(y)):
            if y[i]>=self.fcInitial:                                    #当应力大于初始塑性应力
                if len(xc)==1:                                          #当初始屈服
                    dx = x1[i]-y[i]*self.fc/self.elas/self.ec0          #前期线性段引起应力应变曲线左移量
                    xc.append(x1[i]-dx)                                 #记录为初始屈服点应变
                    yc.append(y[i])                                     #记录为初始屈服点应力
                elif len(xc)>1 and y[i]-yc[-1]>=self.dfct and 1-y[i]>self.dfct/2.0:   #当后屈服点大于前屈服点self.dfc，且小于1-self.dfc/2.0时
                    xc.append(x1[i]-dx)                                 #添加记录点应变
                    yc.append(y[i])                                     #添加记录点应力
                    lineElas = (yc[-1]-yc[-2])/(xc[-1]-xc[-2])          #前后连个点计算切线模量
                    if lineElas>checkElas_c:                            #切线模量需要比初始模量小
                        getWarningReply(message="Error : Data Error!",buttons=(CANCEL,))
                        return 0
        xc.append(1-dx)                                                 #记录峰值应变
        yc.append(1)                                                    #记录峰值应力
        x2 = 1+np.array(range(1,int(((self.ecp)/self.ec0-1+dx)*1000)))/1000.
        y = x2/(self.alpha_d*(x2-1)**2+x2)
        for i in xrange(len(y)):
            if yc[-1]-y[i]>=self.dfct:
                xc.append(x2[i]-dx)
                yc.append(y[i])
            if y[i]<0.01:break
        #判定尾部数据和末尾数据间隔是否大于0.1
        if x2[i]-dx>=xc[-1]+0.1:        #如果大于0.1，则尾部增加数据
            xc.append(x2[i]-dx)
            yc.append(y[i])
        else:                           #否则尾部替换为截断应变
            xc[-1] = x2[i]-dx
            yc[-1] = y[i]
        self.xc = [i*self.ec0 for i in xc]                                          #总压应变
        self.yc = [i*self.fc for i in yc]                                           #总压应力
        self.inpc = [self.xc[i]-self.yc[i]/self.elas for i in range(len(self.xc))]  #非弹性压应变
        self.inpc[0]=0                  #初始非弹性应变为0
        self.inpc[1]=0                  #初始屈服点对应的非弹性应变也为0
        #========================拉应力-应变曲线====================
        xt = [0]
        yt = [0]
        dx = 0                          #前期线性段引起应力应变曲线左移量
        y = 1.2*x1-0.2*x1**6            #归一化应力
        for i in xrange(len(y)):
            if y[i]>=self.ftInitial:
                if len(xt)==1:
                    dx = x1[i]-y[i]*self.ft/self.elas/self.et0
                    xt.append(x1[i]-dx)
                    yt.append(y[i])
                elif len(xt)>1 and y[i]-yt[-1]>=self.dfct and 1-y[i]>self.dfct/2.0:
                    xt.append(x1[i]-dx)
                    yt.append(y[i])
        xt.append(1-dx)
        yt.append(1)
        x2 = 1+np.array(range(1,int(((self.etp)/self.et0-1+dx)*1000)))/1000.
        y = x2/(1.25*np.power((x2-1.0),1.7)+x2)
        for i in xrange(len(y)):
            if yt[-1]-y[i]>=self.dfct:
                xt.append(x2[i]-dx)
                yt.append(y[i])
            if y[i]<0.01:break
        #判定尾部数据和末尾数据间隔是否大于0.1
        if x2[-1]-dx>=xt[-1]+0.1:   #如果大于0.1，则尾部增加数据
            xt.append(x2[i]-dx)
            yt.append(y[i])
        else:                       #否则尾部替换为截断应变
            xt[-1] = x2[i]-dx
            yt[-1] = y[i]
        self.xt = [i*self.et0 for i in xt]                                              #总拉应变
        self.yt = [i*self.ft for i in yt]                                               #总拉应力
        self.inpt = [self.xt[i]-self.yt[i]/self.elas for i in range(len(self.xt))]      #非弹性拉应变
        self.inpt[0]=0              #初始非弹性应变为0
        self.inpt[1]=0              #初始屈服点对应的非弹性应变也为0
        print "cal_Plastic_Brick"
        return 1
    
        
    def cal_DamageEnergy(self):
        #损伤理论，积分形式
        #压缩损伤
        self.dc = [(self.elas*self.xc[i]**2/2.0-(self.cal_dataIntergal(self.xc[:i+1],self.yc[:i+1])))/(
            self.elas*self.xc[i]**2/2.0)  for i in range(1,len(self.xc))]
        #拉伸损伤
        self.dt = [(self.elas*self.xt[i]**2/2.0-(self.cal_dataIntergal(self.xt[:i+1],self.yt[:i+1])))/(
            self.elas*self.xt[i]**2/2.0)  for i in range(1,len(self.xt))]
        #初始损伤需要为0
        self.dc[0] = 0
        self.dt[0] = 0
        #记录压缩残余应变
        self.spc = [0]
        self.spc += [self.xc[i+1]-self.yc[i+1]/(1-self.dc[i])/self.elas for i in range(len(self.dc))]
        self.spc[1] = 0
        #记录拉伸残余应变
        self.spt = [0]
        self.spt += [self.xt[i+1]-self.yt[i+1]/(1-self.dt[i])/self.elas for i in range(len(self.dt))]
        self.spt[1] = 0
    
    def cal_dataIntergal(self,x,y):
        #对数据进行积分求解
        dataSum = 0
        for i in range(len(x)-1):
            dataSum += (x[i+1]-x[i])*(y[i+1]+y[i])/2.0
        return dataSum
    
    def createAbaqusMat(self):
        #在Abaqus中创建材料
        concretemat=mdb.models[self.modelName].Material(name=self.matName)
        concretemat.Elastic(table=((self.elas, self.pos), ))
        concretemat.ConcreteDamagedPlasticity(table=((self.p1, self.p2, self.p3, self.p4, self.p5),))
        #CDP损伤定义
        input_dc=[]
        input_dt=[]
        for i in range(len(self.dc)):
            if self.dc[i]>self.maxDamage:
                break
            input_dc.append([self.dc[i],self.inpc[i+1]])
        for i in range(len(self.dt)):
            if self.dt[i]>self.maxDamage:
                break
            input_dt.append([self.dt[i],self.inpt[i+1]])
        #CDP塑性定义
        table_c = []
        for i in range(len(self.yc)):
            table_c.append([self.yc[i],self.inpc[i]])
        table_t = []
        for i in range(len(self.yt)):
            table_t.append([self.yt[i],self.inpt[i]])
        concretemat.concreteDamagedPlasticity.ConcreteCompressionHardening(table=table_c[1:])
        concretemat.concreteDamagedPlasticity.ConcreteTensionStiffening(table=table_t[1:])
        concretemat.concreteDamagedPlasticity.ConcreteCompressionDamage(table=input_dc,tensionRecovery=self.wt)
        concretemat.concreteDamagedPlasticity.ConcreteTensionDamage(table=input_dt,compressionRecovery=self.wc)
        concretemat.Density(table=((self.den, ),))


#从xuzhouCDP_Brick类中继承绝大多数的功能函数和属性，仅需重新定义参数输入和塑性计算公式
class xuzhouCDP_mortar(xuzhouCDP_Brick):
    def createMaterial(self,matName,den,elas,pos,fc,ft,p1,p2,p3,p4,p5,wt,wc,
    fcInitial,dfct,maxDamage,ec0,ecp,modelName):
        self.matName = matName          #材料名称
        self.den = den                  #密度
        self.elas = elas                #弹性模量
        self.pos = pos                  #泊松比
        self.fc = fc                    #极限压缩应力
        self.ft = ft                    #极限拉伸应力
        self.fcInitial = fcInitial      #初始压缩屈服应力
        self.ec0 = ec0                  #极限压应力对应的应变
        self.ecp = ecp                  #取值最大压应变
        self.et0 = 0.0001               #极限拉应力对应的应变
        self.etp = 0.001                #取值最大拉应变
        self.p1 = p1                    #膨胀角
        self.p2 = p2                    #偏心系数
        self.p3 = p3                    #双向/单向抗压强度比
        self.p4 = p4                    #不变应力比
        self.p5 = p5                    #粘聚系数
        self.wt = wt                    #拉伸恢复系数
        self.wc = wc                    #压缩恢复系数
        self.dfct = dfct                #用于控制输出，当前后两个数据变化大于该数据时设置输出
        self.maxDamage = maxDamage
        self.modelName = modelName      #模型名称
        
        #【计算塑性】
        if self.cal_Plastic()==0:
            return 0
        
        #【计算损伤】
        self.cal_DamageEnergy()
        
        #【数据截断】
        self.yc = [round(i,5) for i in self.yc]
        self.yt = [round(i,5) for i in self.yt]
        self.dc = [round(i,8) for i in self.dc]
        self.dt = [round(i,8) for i in self.dt]
        self.inpc = [round(i,8) for i in self.inpc]
        self.inpt = [round(i,8) for i in self.inpt]
        self.elas=self.elas
        #【生成材料】
        self.createAbaqusMat()
         #修改拉伸砂浆数据
        concretemat=mdb.models[self.modelName].materials[self.matName]
        concretemat.concreteDamagedPlasticity.concreteTensionStiffening.setValues(table=[[self.ft,0],[0.1*self.ft,self.etp]])
        concretemat.concreteDamagedPlasticity.concreteTensionDamage.setValues(table=[[0,0],[0.9,self.etp]])
        print "Finish"
        return 1
    
    def cal_Plastic(self):
        #计算砂浆的塑性
        #========================压应力-应变曲线========================
        xc = [0]
        yc = [0]
        dx = 0                                                          #前期线性段引起应力应变曲线左移量
        x1 = np.array(range(1001))/1000.0                               #归一化应变
        y = x1/(0.3*x1**2+0.4*x1+0.3)                                   #归一化应力
        checkElas_c = self.ec0/(self.fc/self.elas)                      #初始模量的归一化变换
        for i in range(len(y)):
            if y[i]>=self.fcInitial:                                    #当应力大于初始塑性应力
                if len(xc)==1:                                          #当初始屈服
                    dx = x1[i]-y[i]*self.fc/self.elas/self.ec0          #前期线性段引起应力应变曲线左移量
                    xc.append(x1[i]-dx)                                 #记录为初始屈服点应变
                    yc.append(y[i])                                     #记录为初始屈服点应力
                elif len(xc)>1 and y[i]-yc[-1]>=self.dfct and 1-y[i]>self.dfct/2.0:   #当后屈服点大于前屈服点self.dfc，且小于1-self.dfc/2.0时
                    xc.append(x1[i]-dx)                                 #添加记录点应变
                    yc.append(y[i])                                     #添加记录点应力
                    lineElas = (yc[-1]-yc[-2])/(xc[-1]-xc[-2])          #前后连个点计算切线模量
                    if lineElas>checkElas_c:                            #切线模量需要比初始模量小
                        getWarningReply(message="Error : Data Error!",buttons=(CANCEL,))
                        return 0
        xc.append(1-dx)                                                 #记录峰值应变
        yc.append(1)                                                    #记录峰值应力
        x2 = 1+np.array(range(1,int(((self.ecp)/self.ec0-1+dx)*1000)))/1000.
        y = 1.1-0.1*x2
        for i in xrange(len(y)):
            if yc[-1]-y[i]>=self.dfct:
                xc.append(x2[i]-dx)
                yc.append(y[i])
            if y[i]<0.01:break
        #判定尾部数据和末尾数据间隔是否大于0.1
        if x2[i]-dx>=xc[-1]+0.1:    #如果大于0.1，则尾部增加数据
            xc.append(x2[i]-dx)
            yc.append(y[i])
        else:                       #否则尾部替换为截断应变
            xc[-1] = x2[i]-dx
            yc[-1] = y[i]
        self.xc = [i*self.ec0 for i in xc]                                          #总压应变
        self.yc = [i*self.fc for i in yc]                                           #总压应力
        self.inpc = [self.xc[i]-self.yc[i]/self.elas for i in range(len(self.xc))]  #非弹性压应变
        self.inpc[0]=0              #初始非弹性应变为0
        self.inpc[1]=0              #初始屈服点对应的非弹性应变也为0
        #========================拉应力-应变曲线====================
        self.xt = [0,self.ft/self.elas,22*self.ft/self.elas]                        #总拉应变
        self.yt = [0,self.ft,0.1*self.ft]                                           #总拉应力
        self.inpt = [self.xt[i]-self.yt[i]/self.elas for i in range(len(self.xt))]  #非弹性拉应变
        self.inpt[0]=0              #初始非弹性应变为0
        self.inpt[1]=0              #初始屈服点对应的非弹性应变也为0
        print "cal_Plastic_mortar"
        return 1


#判定inpoutput文件夹是否存在
if not os.path.exists(outputDir):
    os.makedirs(outputDir) 

os.chdir(outputDir)
#打开cae文件
openMdb(pathName=caeFullPath)                           #打开cae文件
if len(mdb.models.keys())==1:
    for fci,fc in enumerate(fc_Brick_Lime):
        fc_brick = fc[0]                                #砖抗压强度
        ft_brick = fc[0]/10.0                           #砖抗拉强度
        fc_mortar = fc[1]                               #砂浆抗压强度
        ft_mortar = fc[1]/10.0                          #砂浆抗拉强度
        elas_mortar = fc[2]                             #砂浆弹性模量
        ec0_brick = fc[3]                               #砖极限压应力对应的应变
        ecp_brick = fc[4]                               #砖最大压应变
        ec0_mortar = fc[5]                              #砂浆极限压应力对应的应变
        ecp_mortar = fc[6]                              #砂浆最大压应变

        modelName = mdb.models.keys()[0]
        #创建材料
        brickMain = xuzhouCDP_Brick()
        out1 = brickMain.createMaterial(matName_brick,den_brick,elas_brick,pos_brick,fc_brick,
            ft_brick,p1_brick,p2_brick,p3_brick,p4_brick,p5_brick,wt_brick,wc_brick,
            fcInitial_brick,ftInitial_brick,dfct_brick,maxDamage_brick,alpha_a_brick,
            alpha_d_brick,ec0_brick,ecp_brick,et0_brick,etp_brick,modelName)
        mortarMain = xuzhouCDP_mortar()
        out2 = mortarMain.createMaterial(matName_mortar,den_mortar,elas_mortar,pos_mortar,
            fc_mortar,ft_mortar,p1_mortar,p2_mortar,p3_mortar,p4_mortar,p5_mortar,wt_mortar,wc_mortar,
            fcInitial_mortar,dfct_mortar,maxDamage_mortar,ec0_mortar,ecp_mortar,modelName)
        #创建任务并输出inp
        jobName = "%d-%s_B%.2f_L%.5f"%(fci,modelName,fc_brick/1000000.0,ft_mortar/1000000.0) #任务名称为cae名称+model名称
        jobName = jobName.replace(".","_")          #小数点换下划线
        job = mdb.Job(name=jobName, model=modelName, description='', 
            type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None, 
            memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
            explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
            modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
            scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1, 
            numDomains=1, numGPUs=0)                #创建任务
        job.writeInput(consistencyChecking=OFF)     #输出inp文件
else:
    getWarningReply(message="Only one model is allowed in CAE file!",buttons=(CANCEL,))




# Function: 同时提交多个inp文件；增加时间输出
# 运行方式：1)运行同目录下的calMultiINP.bat文件，并将这两个文件和inp文件放到相同目录下
#===================================================================================================
import time,os
# 开始时间
starttime = time.time()
# 目录文件:根据inp文件路径进行修改
root_dir = "D:/temp/Wall-test-xuzhou/W1/outputInp/"
timeFile = root_dir+"time.csv"
cpuNum = 3      #每个任务占用的cpu数量


# 根据 log 文件判断 job 运行的状态
def statusOfRunningJob(running_job):
    # 读取这个正在运行的任务的日志文件
    print root_dir+ running_job + ".log"
    file = open(root_dir+ running_job + ".log","r")
    content = file.read()
    file.close()
    # 如果完成返回 "finish",失败返回 "error",正在运行返回 "running"
    if "COMPLETED" in content:
        return "finish"
    if "exited with" in content:
        return "error"
    return "running"


#先构造任务列表,失败列表
#任务列表为当前目录下的所有inp文件
job_name_list = []
for item in os.listdir(root_dir):
    ext=os.path.splitext(item)[1]               #后缀
    if ext == ".inp":                           #判断是否为inp文件
        odbName = item.replace(".inp",".odb")
        if not os.path.isfile(odbName):         #判断是否存在同名odb，如果存在则表示已经计算
            job_name_list.append(item[:-4])

print job_name_list
error_job_list = []     #失败列表


# 再构造运行队列
running_job_queue = []
# 设置队列最大长度
max_len = 4

# 循环遍历任务列表
for job_name in job_name_list:
    # 先检查运行队列的状态,如果队列已满,进入循环
    while len(running_job_queue) == max_len:
        #   循环遍历正在运行的任务
        for running_job in running_job_queue:
            # 查看该任务运行状态
            tag = statusOfRunningJob(running_job)
            # 如果任务不运行了,处于完成或者错误状态
            # 就从 running 列表中剔除它
            if tag != "running":
                running_job_queue.remove(running_job)
                # 如果任务失败了,我们还要用 error_job_list 来记录它
                if tag == "error":
                    error_job_list.append(running_job)
        # 每一次循环需要停 5s ,不然 cpu 全用在循环上了.
        time.sleep(5)
        print job_name
    # 直到队列不非满时,
    # 1 通过 abaqus 接口提交任务
    print "abaqus job=%s cpus=%d"%(job_name,cpuNum)
    os.popen("abaqus job=%s cpus=%d"%(job_name,cpuNum))
    # 2 把任务名加入到running_job_queue中
    running_job_queue.append(job_name)

#最后几个任务需要监控
while len(running_job_queue) >0:
    for running_job in running_job_queue:
        tag = statusOfRunningJob(running_job)
        if tag != "running":
            running_job_queue.remove(running_job)
            if tag == "error":
                error_job_list.append(running_job)
    time.sleep(5)
    print len(running_job_queue)


# 记录错误任务名到 error.log 文件中
file = open("error.log", "w")
file.write("\n".join(error_job_list))
file.close()
# 结束时间
endtime = time.time()
total_time = endtime - starttime 
#print "总程序运行时间为{}".format(total_time)
with open(timeFile, 'a+') as fp:
    fp.write(str(total_time)+"\n")


# Function: 处理滞回曲线数据获得骨架曲线,增加提取结构文件目录，并避免重复提取
#================================================================================
from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
import displayGroupOdbToolset as dgo
import customKernel
#获得指定目录下的odb文件名称列表
odbFilePath = "D:/temp/Wall-test-xuzhou/W1/outputInp/"            #需要修改的目录名称，名称后面注意有个反斜杠
outPutPath = "D:/temp/Wall-test-xuzhou/W1/outputInp/result/"      #指定输出目录
if not os.path.exists(outPutPath):                                #判断工作目录是否存在
    os.makedirs(outPutPath)                                       #递归创建目录

odbNames = []
for item in os.listdir(odbFilePath):                              #遍历目录下的所有文件
    ext=os.path.splitext(item)[1]                                 #后缀
    if ext == ".odb":                                             #判断是否为odb文件
        odbNames.append(odbFilePath+item)                         #如果是，则存储文件路径和名称

print odbNames


outPutMaxDataFile = outPutPath+"maxData.csv"
if os.path.isfile(outPutMaxDataFile):                             #判断路径是否为文件
    f = open(outPutMaxDataFile,"r")
    fStr = f.readlines()
    f.close()
    fmaxData = open(outPutMaxDataFile,"w")
    #fmaxData.write(fStr[0])
    for line in fStr[1:]:
        odbName = line.split(",")[0]
        fmaxData.write(line)
        if odbName in odbNames:
            odbNames.remove(odbName)
else:
    fmaxData = open(outPutMaxDataFile,"w")
    fmaxData.write("odbName,maxForce,maxDisplacement")

fmaxData.close()
print len(odbNames)



def getRockDatas (curOdbName):
    #获得骨架曲线
    curOdbName
    odb = session.openOdb(curOdbName)                                                     #odb对象
    rpSet = "rp1"                                                                         #参考点集合名称
    stepName = "Step-3"                                                                   #提取数据的分析步
    times = []                                                                            #存储时间
    u1 = []                                                                               #存储位移
    rf1 = []                                                                              #存储反力
    times_u1 = []                                                                         #存储时间位移数据
    times_rf1 = []                                                                        #存储时间反力数据
    u1_rf1 = []                                                                           #存储滞回曲线数据
    region = odb.rootAssembly.nodeSets[rpSet.upper()]                                     #参考点集合
    csvName = outPutPath + os.path.split(curOdbName)[-1].replace(".odb",".csv")           #存储滞回数据文件
    rockFileName = outPutPath + os.path.split(curOdbName)[-1].replace(".odb","-rock.csv") #存储骨架数据文件
    rockPngName = curOdbName.replace(".odb","-rock.png")        #存储png图片数据
    #从odb中读取数据
    stepFrames = odb.steps[stepName].frames                     #提取的时间帧为最后一个分析步
    for frame in stepFrames:                                    #遍历时间
        frameTime = frame.frameValue                            #帧时间
        frameU1 = frame.fieldOutputs['U'].getSubset(
            region = region).values[0].data[0]                  #u1数值
        frameRF1 = frame.fieldOutputs['RF'].getSubset(
            region = region).values[0].data[0]                  #rf1数值
        #存储数据
        times.append(frameTime)
        u1.append(frameU1)
        rf1.append(frameRF1)
        times_u1.append([frameTime,frameU1])
        times_rf1.append([frameTime,frameRF1])
        u1_rf1.append([frameU1,frameRF1])
    
    #存储到滞回数据文件
    f = open(csvName,"w")                                       #打开文件
    f.write("time,u1,rf1")                                      #写入表头
    for i in range(len(times)):                                 #循环写入
        f.write("\n%f,%f,%f"%(times[i],u1[i],rf1[i]))           #写入数据
    f.close()                                                   #关闭文件
    
    
    #绘制曲线
    curve1 = session.XYData(name="times_u1",data=times_u1)      #times_u1曲线
    curve2 = session.XYData(name="times_rf1",data=times_rf1)    #times_rf1曲线
    c1 = session.Curve(xyData=curve1)
    c2 = session.Curve(xyData=curve2)
    xQuantity = visualization.QuantityType(type=DISPLACEMENT)
    yQuantity = visualization.QuantityType(type=FORCE)
    curve3 = session.XYData(name="u1_rf1",data=u1_rf1,
        xValuesLabel="Displacement",yValuesLabel="Force",
        axis1QuantityType=xQuantity,axis2QuantityType=yQuantity)#滞回曲线
    c3 = session.Curve(xyData=curve3)
    c3.lineStyle.setValues(style=DOTTED,color='#0000FF')        #设置滞回曲线线型和颜色
    if 'XYPlot-1' in session.xyPlots.keys():
        xyp = session.xyPlots['XYPlot-1']
    else:
        xyp = session.XYPlot(name='XYPlot-1')
    chart = xyp.charts[xyp.charts.keys()[0]]


    #数据处理获得反力极值点位置
    indexs = [0]
    Points = [[0,0]]
    rf1.append(0)
    for i in range(1,len(times)-1):
        d1 = rf1[i-1]
        d2 = rf1[i]
        d3 = rf1[i+1]
        if(d2-d1)*(d2-d3)>0:                                    #通过前后两个点的差值乘积判断是否为极值点
            indexs.append(i)
            Points.append([times[i],rf1[i]])
    Points.append([times[-1],0])

    curve2Max = session.XYData(name="times_rfMax",data=Points)
    c2Max = session.Curve(xyData=curve2Max)
    c2Max.symbolStyle.setValues(show=True,size=3)
    c2Max.lineStyle.setValues(show=False)

    #剔除小波峰
    newPoints = Points[:]
    limit = abs(Points[1][1]*0.9)   #小波峰的筛选条件设置为第一个波峰的0.9倍
    mark = 1                        #用于记录是否存在小波峰
    while mark:
        mark = 0
        minLen = limit
        minLenId = []
        #遍历查找相邻两个极值点的距离，如果为最小，且小于limit，则需要剔除
        for i in range(0,len(newPoints)-1):
            d1= newPoints[i][1]
            d2= newPoints[i+1][1]
            if abs(d1-d2)<minLen:
                mark=1;
                minLen = abs(d1-d2)
                minLenId = [i,i+1]
        if mark:
            newPoints.remove(newPoints[minLenId[1]])
            newPoints.remove(newPoints[minLenId[0]])
            indexs.remove(indexs[minLenId[1]])
            indexs.remove(indexs[minLenId[0]])
    
    curve2Max2 = session.XYData(name="times_rfMax2",data=newPoints)
    c2Max2 = session.Curve(xyData=curve2Max2)
    c2Max2.symbolStyle.setValues(show=True,size=2)
    c2Max2.lineStyle.setValues(show=False)
    
    #骨架曲线
    rockData1 = [[0,0]]
    rockData2 = [[0,0]]
    #极大值和极小值是间隔出现，因此交替存储
    for i in range(len(indexs)):
        if i%2:
            rockData1.append([u1[indexs[i]],rf1[indexs[i]]])
        else:
            rockData2.append([u1[indexs[i]],rf1[indexs[i]]])
    rockData2.reverse()                                         #将极小值倒叙
    rockData = rockData2+rockData1                              #合并为整个连续的骨架曲线
    rock = session.XYData(name="rock",data=rockData,
        xValuesLabel="Displacement",yValuesLabel="Force",
        axis1QuantityType=xQuantity,axis2QuantityType=yQuantity)
    rockCurve = session.Curve(xyData=rock)
    rockCurve.symbolStyle.setValues(show=True,size=2,marker=FILLED_SQUARE,color='#FF0000')
    rockCurve.lineStyle.setValues(color='#FF0000')
    chart.setValues(curvesToPlot=(c3), )
    curView = session.viewports[session.currentViewportName]
    curView.setValues(displayedObject=xyp)
    
    #输出数据
    f = open(rockFileName,"w")
    f.write("u1,rf1\n")
    maxForce = 0
    maxDisplacement = 0
    for i in rockData:
        f.write("\n%f,%f"%(i[0],i[1]))
        if abs(i[1])>maxForce:
            maxForce=abs(i[1])
            maxDisplacement = i[0]
    f.close()
    
    #输出图片
    session.printOptions.setValues(vpDecorations=OFF)
    session.printToFile(fileName=rockPngName, format=PNG, canvasObjects=(curView, ))
    odb.close()
    
    #输出极大值
    f = open(outPutMaxDataFile,"r")
    fStr = f.read()
    f.close()
    fmaxData = open(outPutMaxDataFile,"w")
    fmaxData.write(fStr)
    fmaxData.write("\n%s,%f,%f"%(curOdbName,maxForce,maxDisplacement))
    fmaxData.close()

#循环提取
for curOdbName in odbNames:
    getRockDatas(curOdbName)