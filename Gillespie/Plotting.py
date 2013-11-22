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