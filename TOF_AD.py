
from uniqed.runners.tof_run import detect_outlier
import matplotlib.pyplot as plt
import pandas as pd
import string



#print(data_df[['connectance']])
#Detect Outlier method adapted from uniqed
def TOF_ADector(data_list:list, s_index:int,  variable_name: string, non_s:bool=True):
    input_data=data_list[s_index][[variable_name]]
    #pd.set_option('display.max_rows', None)

    if non_s==False:

        preprocessed_data=input_data.diff()
        preprocessed_data.iloc[0]=0
        res_df = detect_outlier(preprocessed_data, cutoff_n=10)
    #print(res_df)
    else:
        res_df = detect_outlier(input_data, cutoff_n=10)

    fig, axs = plt.subplots(2, 1, sharex=True)
    fig.suptitle(f'TOF anomaly detection of {variable_name} (sample {s_index+1})')

    if non_s==False:

        axs[0].plot(res_df[variable_name], color='tab:blue', label='diff time series')
        axs[0].set_ylabel(f"diff_{variable_name}")
    else:
        axs[0].plot(res_df[variable_name], color='tab:blue', label='time series')
        axs[0].set_ylabel(variable_name)
    #axs[0].plot(res_df['connectance'].loc[data_df.query("is_anomaly==1").index.values],
    #         color='tab:green', label='anomaly')
    axs[0].plot(res_df.query("TOF==1")[variable_name], lw=0, marker='o', ms=0.5,
             color='tab:orange', label='TOF detections')

    axs[0].legend(loc='upper right', framealpha=0.6)


    axs[1].plot(res_df['TOF_score'], color='k', label='TOF score')
    #axs[1].plot(res_df['TOF_score'].loc[data_df.query("is_anomaly==1").index.values],
    #         color='tab:green', label='anomaly')
    axs[1].plot(res_df.query("TOF==1")['TOF_score'], lw=0, marker='o', ms=0.5,
             color='tab:orange', label='TOF')
    axs[1].set_ylabel('TOF score')
    axs[1].set_xlabel('t')
    axs[1].legend(['TOF score','TOF detections'],
                  loc='upper right',
                  framealpha=0.6)

    axs[1].set_xlim(0, 600)
    axs[0].grid(True)
    axs[1].grid(True)

    if non_s==True:
    #fig.tight_layout(rect=[0, 0, 1, 1], pad=1, h_pad=0, w_pad=0)
        fig.savefig(f"{variable_name}_sample{s_index+1}_run.png")

    else:
        fig.savefig(f"stationary_{variable_name}_sample{s_index+1}_run.png")




    print(res_df.query("TOF==1")['TOF_score'])
    plt.show()

#for decomposed datasets
def TOF_ADector_Decomp(data_list:list, s_index:int, variable_name: string):
    input_data=data_list[s_index][['remainder']]


    res_df = detect_outlier(input_data, cutoff_n=10)

    fig, axs = plt.subplots(3, 1, sharex=True)
    fig.suptitle(f'TOF anomaly detection of {variable_name} (sample {s_index+1})')

    axs[0].plot(data_list[s_index][['observed']], color='tab:blue', label='time series')

    axs[1].plot(res_df['remainder'], color='tab:red', label='remainder of time series')

    #axs[0].plot(res_df['connectance'].loc[data_df.query("is_anomaly==1").index.values],
    #         color='tab:green', label='anomaly')
    axs[1].plot(res_df.query("TOF==1")['remainder'], lw=0, marker='o', ms=2.5,
             color='tab:orange', label='TOF detections')
    axs[0].set_ylabel(variable_name)
    axs[0].legend(loc='upper right', framealpha=0.6, fontsize="x-small")

    axs[1].set_ylabel('remainder')
    axs[1].legend(loc='upper right', framealpha=0.6, fontsize="x-small")


    axs[2].plot(res_df['TOF_score'], color='k', label='TOF score')
    #axs[1].plot(res_df['TOF_score'].loc[data_df.query("is_anomaly==1").index.values],
    #         color='tab:green', label='anomaly')
    axs[2].plot(res_df.query("TOF==1")['TOF_score'], lw=0, marker='o', ms=2,
             color='tab:orange', label='TOF')
    axs[2].set_ylabel('TOF score')
    axs[2].set_xlabel('t')
    axs[2].legend(['TOF score','TOF detections'],
                  loc='upper right',
                  framealpha=0.6,
                  fontsize='x-small')

    axs[2].set_xlim(0, 600)
    axs[0].grid(True)
    axs[1].grid(True)
    axs[2].grid(True)



    fig.savefig(f"Decomposed_{variable_name}_sample{s_index+1}_run.png")




    print(res_df.query("TOF==1")['TOF_score'])
    plt.show()



def main():
    #Data Preparation
    df_list=[]
    Conn_list_rem=[]
    modu_list_rem=[]
    cc_list_rem=[]
    nest_list_rem=[]
    compart_list_rem=[]


    for i in range(1, 11):
        df_list.append(pd.read_csv(f"indicators_sample{i}.csv"))
        Conn_list_rem.append(pd.read_csv(f"Connectance_sample{i}_decomposed.csv"))
        modu_list_rem.append(pd.read_csv(f"Modu_sample{i}_decomposed.csv"))
        cc_list_rem.append(pd.read_csv(f"CC_sample{i}_decomposed.csv"))
        nest_list_rem.append(pd.read_csv(f"Nest_sample{i}_decomposed.csv"))
        compart_list_rem.append(pd.read_csv(f"#Compartment_sample{i}_decomposed.csv"))



    #print(df_list[2]['connectance'][200:300])
    #Detection
    #TOF_ADector(df_list, 0, 'connectance', non_s=False)

    TOF_ADector_Decomp(Conn_list_rem, 5, "connectance")
    TOF_ADector_Decomp(modu_list_rem, 5, "modularity")
    TOF_ADector_Decomp(cc_list_rem, 5, 'clustering coefficient')
    #TOF_ADector_Decomp(nest_list_rem, 0, "nestedness")
    #TOF_ADector_Decomp(compart_list_rem, 4, "# of Compartment")


if __name__=="__main__":
    main()
