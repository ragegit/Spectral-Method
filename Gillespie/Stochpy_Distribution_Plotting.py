#use this to extract data:
data_m = smod.data_stochsim.getSimData('m')
data_n = smod.data_stochsim.getSimData('n')

distributions = smod.data_stochsim.species_distributions


def GetDataDistributions(L_sim_output,L_names):
    """
    Get distributions, means, standard deviations, and the (raw) moments
      
    Input:
    - *L_sim_output* (list)
    - *L_names* (list)
    """
    n_names = len(L_names)
    L_distributions = [{} for i in xrange(n_names)]
    endtime = L_sim_output[-1][0]
    for timestep in xrange(len(L_sim_output)-1):
        for i in xrange(n_names):
            try:
                L_distributions[i][L_sim_output[timestep][i+1]] += L_sim_output[timestep+1][0] - L_sim_output[timestep][0]
            except KeyError:
                L_distributions[i][L_sim_output[timestep][i+1]] = L_sim_output[timestep+1][0] - L_sim_output[timestep][0]

    L_probability_mass = []
    D_means = {}
    D_stds = {}
    D_moments = {}
    for i in xrange(n_names):
        L_x_i = np.array(sorted(L_distributions[i]))
        L_y_i = np.array([L_distributions[i][x_i] for x_i in sorted(L_distributions[i])])/float(endtime) # probability = dt/T
        mean = (L_x_i*L_y_i).sum()
        mean_sq = (L_x_i**2*L_y_i).sum()
        var = mean_sq - mean*mean
        std = var**0.5
        L_probability_mass.append([L_x_i,L_y_i])
        id_name = L_names[i]

        D_means[id_name] = mean
        D_stds[id_name] = std

        D_moments[id_name] = {}
        D_moments[id_name]['1'] = mean
        D_moments[id_name]['2'] = mean_sq
        D_moments[id_name]['3'] = (L_x_i**3*L_y_i).sum()     
        D_moments[id_name]['4'] = (L_x_i**4*L_y_i).sum()
    del L_sim_output,L_distributions
    return (L_probability_mass,D_means,D_stds,D_moments)


      def AverageDistributions(self,L_means,L_stds,nstd,datatype2plot,L_labels,linestyle,marker_,colors,title,xlabel,ylabel,IsLegend):
      """      
      Plots the average and standard deviation
      
      Input:
       - *L_means* (nested list)
       - *L_stds* (nested list)
       - *nstd* (float)
       - *L_labels* (list)
       - *linestyle* (string)
       - *title* (string)
       - *xlabel* (string)
       - *ylabel* (string)
       - *IsLegend* (boolean)
      """
      assert nstd > 0, "Error: The number of STDs must be a value larger than zero"
      plt.figure(self.plotnum)
      datatype2plot_indices = [L_labels.index(datatype) for datatype in datatype2plot]
      j=0
      for i in datatype2plot_indices:
          if colors == None:
              if j >= len(self.colors):
                  j=0
          elif colors:
              if j >= len(colors):
                  j=0
          if colors == None:
              plt.errorbar(L_means[i][0],L_means[i][1],yerr = nstd * np.array(L_stds[i][1]),color = self.colors[j],ls = linestyle,marker = marker_,label = L_labels[i]) # plot with y-axis error bars
          else:
              if clr.is_color_like(colors[j]):
                  plt.errorbar(L_means[i][0],L_means[i][1],yerr = nstd*np.array(L_stds[i][1]),color = colors[j],ls = linestyle,marker = marker_,label = L_labels[i])
              else:
                  print "Warning: '%s' is not recognized as a valid color code" % colors[j]
                  plt.errorbar(L_means[i][0],L_means[i][1],yerr = nstd * np.array(L_stds[i][1]),color = self.colors[j],ls = linestyle,marker = marker_,label = L_labels[i])
                  colors = None
          j+=1
      if IsLegend:
          plt.legend(numpoints=1,frameon=True)
      plt.title(title)
      plt.xlabel(xlabel)
      plt.ylabel(ylabel)