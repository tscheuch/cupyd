# Load SWMM model
dire="SWMM_Llanquihue"
file="/Llanquihue_base.inp"
SWMM_path=dire + "/" + file

## Import SWMM objects
with Simulation(SWMM_path) as sim:
    s_list=[]
    su_list=[]
    j_list=[]
    l_list=[]
    s_areas=[]
    for s in s_names_list: #lista de strings de los nombres de las subcuencas en .inp
        s_list.append(Subcatchments(sim)[s]) 
        s_areas.append(Subcatchments(sim)[s].area) #[ha]  
    for su in su_names_list)): #lista de strings de los nombres de las unidades de almacenamiento en .inp
        su_list.append(Nodes(sim)[su])
    for j in j_names_list: #lista de strings de los nombres de los nodos (junctions, outfalls and dividers) en .inp
        j_list.append(Nodes(sim)[j])
    for l in l_names_list: #lista de strings de los nombres de los conductos (links) en .inp
        l_list.append(Links(sim)[l]) 

## Create empty DataFrames for SWMM objects Time Series Results: 
s_SWMM_ts=[]
su_SWMM_ts=[]
j_SWMM_ts=[]
l_SWMM_ts=[]
for i in s_names_list:
    a= pd.DataFrame()
    a["Time"]=""
    a["Precipitation (mm/h)"]=""
    a["Evaporation (mm/d)"]=""
    a["Infiltration (mm/h)"]="" 
    a["Runoff (m3/s)"]=""
    a["Runon (m3/s)"]=""
    a["Cumulative Infiltration (m3)"]=""
    a["Cumulative Evaporation (m3)"]=""
    s_SWMM_ts.append(a)

for i in su_names_list:
    a= pd.DataFrame()
    a["Time"]=""
    a["Depth (m)"]=""
    a["Head (m)"]=""
    a["Flooding (m3/s)"]=""
    a["Lateral inflow (m3/s)"]=""
    a["Total inflow (m3/s)"]=""
    a["Total outflow (m3/s)"]=""
    a["Volume (m3)"]=""
    a["Losses (m3/s)"]=""
    a["Cumulative Exfiltration Loss (m3)"]=""
    a["Cumulative Evaporation Loss (m3)"]=""
    su_SWMM_TS.append(a)

for i in j_names_list:
    a= pd.DataFrame()
    a["Time"]=""
    a["Depth (m)"]=""
    a["Head (m)"]=""
    a["Flooding (m3/s)"]=""
    a["Lateral inflow (m3/s)"]=""
    a["Total inflow (m3/s)"]=""
    a["Total outflow (m3/s)"]=""
    J_SWMM_TS.append(a)
for i in l_names_list:
    a= pd.DataFrame()
    a["Time"]=""
    a["Depth (m)"]=""
    a["Flow (m3/s)"]=""
    l_SWMM_ts.append(a)

# Loop for Temporal/Spatial Data Exchange
TS_gdf=[] # empty matrix for gdf with cells results for each time step
ZB_TS=[] # empty matrix for MODFLOW Zone Budget results. (util para el post proceso)

%%time
%%capture
with Simulation(SWMM_path) as sim:
    
    system_routing = SystemStats(sim)
    #Lists of cummulative infiltration to calculate delta infiltratation for every time step: S and SU 
    inf_s_list_1=np.zeros(len(s_list))
    inf_s_list_2=np.zeros(len(s_list))
    inf_su_list_1=np.zeros(len(su_list))
    inf_su_list_2=np.zeros(len(su_list))

    
    #Lists for the DRN incorporation
    rate_j_list=np.zeros(len(j_list))
    rate_su_list=np.zeros(len(su_list))
  
    
    #time counter for daily infiltration agregation and hourly reports #esto debiera quedar en función del paso de tiempo de MODFLOW y no siempre diario
    step_counter=0
    day_counter=0
    hourly_counter=0
    for step in sim:
        step_counter=step_counter+1

        if step_counter==360: #CHANGE ACCORDING TO SWMM DT
            step_counter=0
            hourly_counter+=1 #CHANGE ACCORDING TO SWMM DT (REPORT)
            
            # TIME SERIES RESULTS: 1HR AS REPORT TIME STEP
            
            #SUBCATCHMENTS TIME SERIES RESULTS
            for i in range(len(s_list)):
                new_row = {'Time': sim.current_time, "Precipitation (mm/h)":s_list[i].rainfall, 
                           "Evaporation (mm/d)":s_list[i].evaporation_loss,"Infiltration (mm/h)":s_list[i].infiltration_loss, 
                           "Runoff (m3/s)":s_list[i].runoff,"Runon (m3/s)":s_list[i].runon,
                           "Cumulative Infiltration (m3)": s_list[i].statistics["infiltration"], "Cumulative Evaporation (m3)": s_list[i].statistics["evaporation"]}
                s_SWMM_ts[i] =s_SWMM_ts[i].append(new_row, ignore_index=True)

            #STORAGE UNITS TIME SERIES RESULTS
            for i in range(len(su_list)):
                new_row = {'Time': sim.current_time, "Depth (m)":su_list[i].depth, 
                           "Head (m)":su_list[i].head, "Flooding (m3/s)":su_list[i].flooding, 
                           "Lateral inflow (m3/s)":su_list[i].lateral_inflow,"Total inflow (m3/s)":su_list[i].total_inflow,
                           "Total outflow (m3/s)":su_list[i].total_outflow, "Volume (m3)":su_list[i].volume, "Losses (m3/s)":su_list[i].losses,
                           "Cumulative Exfiltration Loss (m3)": su_list[i].storage_statistics["exfil_loss"], "Cumulative Evaporation Loss (m3)": su_list[i].storage_statistics["evap_loss"]}
                su_SWMM_TS[i] =su_SWMM_TS[i].append(new_row, ignore_index=True)
            
            #NODES (junctions, dividers and outfalls) TIME SERIES RESULTS
            for i in range(len(j_list)):
                new_row = {'Time': sim.current_time, "Depth (m)":j_list[i].depth, 
                          "Head (m)":j_list[i].head, "Flooding (m3/s)":j_list[i].flooding, 
                          "Lateral inflow (m3/s)":j_list[i].lateral_inflow,"Total inflow (m3/s)":j_list[i].total_inflow,
                          "Total outflow (m3/s)":j_list[i].total_outflow}
                j_SWMM_ts[i] =j_SWMM_ts[i].append(new_row, ignore_index=True) 
                
                
            #CONDUITS TIME SERIES RESULTS
            for i in range(len(l_list)):
                new_row = {'Time': sim.current_time, "Depth (m)":l_list[i].depth, 
                           "Flow (m3/s)":l_list[i].flow}
                l_SWMM_ts[i] =l_SWMM_ts[i].append(new_row, ignore_index=True)
                  
        #DAILY INFILTRATION ON STORAGE UNITS and SUBCATCHMENTS:
        
        if hourly_counter==24: #CHANGE ACCORDING TO MODFLOW DT
            day_counter=day_counter+1
            hourly_counter=0
            
            print(sim.current_time)

            for i in range(len(s_list)):
                #Delta infiltration
                
                inf_s_list_2[i]=(s_list[i].statistics["infiltration"]-inf_s_list_1[i])
                inf_s_list_1[i]=s_list[i].statistics["infiltration"]
                
            for i in range(len(su_list)):
                #Delta infiltration
                inf_su_list_2[i]=(su_list[i].storage_statistics["exfil_loss"]-inf_su_list_1[i])
                inf_su_list_1[i]=su_list[i].storage_statistics["exfil_loss"]

            RCH_s=inf_s_list_2
            RCH_su=inf_su_list_2

            #CHANGE OF UNITS m3/day->m/day:
    
            RCH_s_M=RCH_s/(np.array(s_areas_modflow))
            RCH_su_M=RCH_su/(np.array(su_areas_modflow))
           
            RCH_s_df=pd.DataFrame({"s":s_names_list, "RCH_s":RCH_s_M})
            RCH_su_df=pd.DataFrame({"su":su_names_list, "RCH_su":RCH_su_M})
          

            #Gereferenced RCH: Add to the MODFLOW_gdf_test new columns RCH_s, RCH_su
            
            MODFLOW_gdf_loop=pd.DataFrame()
            
            MODFLOW_gdf_loop=pd.merge(MODFLOW_gdf_SWMM, RCH_s_df, on="s", how="left")
            MODFLOW_gdf_loop=pd.merge(MODFLOW_gdf_loop, RCH_su_df, on="su", how="left")
        
            
            # Sum georeferences RCHs
            
            MODFLOW_gdf_loop["RCH"]= MODFLOW_gdf_loop.RCH_s.fillna(0) + MODFLOW_gdf_loop.RCH_su.fillna(0)
            
            #Create MODFLOW inputs: RCH package
            
            rch_array = np.zeros((nrows,ncols))
            recharge_cells = MODFLOW_gdf_SWMM.index.values
            for cell in recharge_cells:
                row = MODFLOW_gdf_loop.row[cell]
                col = MODFLOW_gdf_loop.column[cell]
                flux = MODFLOW_gdf_loop.RCH[cell]
                rch_array[row - 1][col - 1] = flux
            rch_array[np.isnan(rch_array)] = 0.
        
            rch = flopy.modflow.ModflowRch(ml, nrchop=3,rech=rch_array, ipakcb=53)
            
            ########################################################################
            #Update boundry conditions for lake and river
            
            MODFLOW_gdf_loop["strt"]=""
            strt_df=np.reshape(strt, len(MODFLOW_gdf_loop))
            MODFLOW_gdf_loop["strt"]=strt_df
            for i in range(len(MODFLOW_gdf_loop.strt)):
                if MODFLOW_gdf_loop["river"][i]:
                    MODFLOW_gdf_loop["strt"][i]=boundry_conditions_df.Elev_maullin[sim.current_time]
                if MODFLOW_gdf_loop["lake"][i]:
                    MODFLOW_gdf_loop["strt"][i]=boundry_conditions_df.Elev_llanq[sim.current_time]
            #########################################################################

            strt_array = np.zeros((nrows,ncols))
            strt_cells = MODFLOW_gdf_SWMM.index.values
            for cell in strt_cells:
                row = MODFLOW_gdf_loop.row[cell]
                col = MODFLOW_gdf_loop.column[cell]
                flux = MODFLOW_gdf_loop.strt[cell]
                strt_array[row - 1][col - 1] = flux
            
            #bas = flopy.modflow.ModflowBas(ml, ibound=ibound, strt=strt)
            bas = flopy.modflow.ModflowBas(ml, ibound=ibound, strt=strt_array) #use the head table of the last time step and bc
           
            #Run MODFLOW
            
            ml.write_input()
            ml.run_model(silent=True)
            
            #Read MODFLOW outputs
            fname = os.path.join('LLANQUIHUE_MF_workspace', 'LLANQUIHUE.hds')
            headfile = flopy.utils.HeadFile(fname, model=ml)
            heads = headfile.get_data()
            heads[heads==1.e+30] = np.nan            # fix masked data 
            heads[heads==-999.99] = np.nan

            #Strt next loop
            
            strt = heads[0]
            
            top = ml.dis.top.array
            DTWT = (top-DRN_burn_depth)-heads[0]  
            
            #DRN calculation
            
            delta_H=np.reshape(DTWT, len(MODFLOW_gdf_loop))
            altura=np.reshape(heads[0], len(MODFLOW_gdf_loop))
            for i in range(len(delta_H)):
                if delta_H[i]<0:
                    delta_H[i]=-delta_H[i]
                else:
                     delta_H[i]=0.
            MODFLOW_gdf_loop["Altura"]=altura   
            MODFLOW_gdf_loop["delta_H"]=delta_H
            MODFLOW_gdf_loop["DRN_rate"]=0.
            
            for i in range(len(MODFLOW_gdf_loop["DRN"])):
                if MODFLOW_gdf_loop["DRN"][i]==1:
                    MODFLOW_gdf_loop["DRN_rate"][i]=(delta_H[i])*DRN_C
                if MODFLOW_gdf_loop["ibound"][i]==-1:
                    MODFLOW_gdf_loop["DRN_rate"][i]=0
                
            
            #INFLOW RATES IN SU AND JUNCTIONS
            
            #Inflow rates in SU:
            
            inflow_su=MODFLOW_gdf_loop.groupby("drn_to").sum()["DRN_rate"]
            inflow_su_list=[]
            for su in su_names_list:
                if su in inflow_su.index:
                    inflow_su_list.append(inflow_su[w])
                else:
                    inflow_su_list.append(0.)
    
                    
            #Inflow rate in Junctions
            
            inflow_j=MODFLOW_gdf_loop.groupby("drn_to").sum()["DRN_rate"]
            inflow_j_list=[]
            for j in j_names_list:
                if j in inflow_j.index:
                    inflow_j_list.append(inflow_j[j])
                else:
                    inflow_j_list.append(0.)
                
            #Generated inflow SWMM WSU:
            
            for i in range(len(su_list)):
                rate=inflow_su_list[i]/86400. #m3/s
                su_list[i].generated_inflow(rate)
                
                
            #Generated inflow SWMM Junctions:
            
            for i in range(len(j_list)):
                rate=inflow_j_list[i]/86400 #m3/s
                j_list[i].generated_inflow(rate)
                
            #Save MODFLOW cells information
            TS_gdf.append(MODFLOW_gdf_loop)
            
                
            #Save zonebudget DataFrame
            fname = os.path.join('LLANQUIHUE_MF_workspace', 'LLANQUIHUE.cbc')
            cbb = flopy.utils.CellBudgetFile(fname)
            zb = flopy.utils.ZoneBudget(cbb, zon, aliases=aliases)
            zb_df=zb.get_dataframes()
            ZB_TS.append(zb_df)
            
                
    routing_stats=system_routing.routing_stats
    runoff_stats=system_routing.runoff_stats
            
print("Flow Routing Mass Balance Error:", sim.flow_routing_error)
print("Runoff Mass Balance Error:", sim.runoff_error)

