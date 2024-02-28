# AA_LocationManager.R
# By: J.M. Rodriguez
# Date: 12024-01-19

# Description:
# Manage locations of for the entire set of scripts so that they are all
# referring to the same place.

# --- Script location (working directory) -------------------------------------
if(
  # Verify current working directory
  (file.exists("./Scripts/verifyScripts.txt") &&
       readLines("./Scripts/verifyScripts.txt")[2] == "verifying-12345")){
    
  # If initial wd failed, back up a layer, try again
  }else if((file.exists("../Scripts/verifyScripts.txt") &&
            readLines("../Scripts/verifyScripts.txt")[2] == "verifying-12345")){
    
    setwd("..") # Back up a folder
    
  }else{
    # try getting location of editor. May fail if doing stuff on OSC
    setwd(dirname(rstudioapi::getSourceEditorContext()$path))
    setwd("..")
    
    # Check that editor gave the proper working directory.
    if(!(file.exists("../Scripts/verifyScripts.txt") &&
        readLines("../Scripts/verifyScripts.txt")[2] == "verifying-12345")){
      
      # throw an error if it all ends up failing
      stop("Could not verify location of scripts/working directory...\n",
           "Current working directory set to:\n",
           getwd())
    }
  }

# --- Data location  ----------------------------------------------------------
loc_Gdr <- c(getwd(),
             "G:/Shared drives/NASA_ForestScaling/LifeHistoryStrats/")
i <- 1
repeat{
  if(file.exists(paste0(loc_Gdr[i],"verifyGdr.txt")) &&
     # Verify file "password"
     readLines(paste0(loc_Gdr[i],"verifyGdr.txt"))[2] == "verifying-12345"){
    
    loc_Gdr <- loc_Gdr[i] # Change location of GDrive upon success
    break
    
  }else if(i >= length(loc_Gdr)){ 
    # Ensure that there are more location to check, if not throw error
    stop("'verifyGdr.txt' file not found in any of the provided locations. ",
         "Please ensure that the location of data is set up properly.\n",
         "Locations searched:\n- ", paste(loc_Gdr, collapse = "\n- "))
  }
  
  # if not found, continue to the next location
    i <- i + 1
}

# Notify user where the data is being pulled from
cat("Working directory for scripts set to:\n",
    getwd(),
    "\nLocation of data set to:\n",
    loc_Gdr,
    "\n")

# End