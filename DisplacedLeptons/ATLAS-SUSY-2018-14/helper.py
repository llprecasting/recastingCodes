#!/usr/bin/env python3
import numpy as np
from pathlib import Path
from scipy.interpolate import NearestNDInterpolator
from numpy import ndarray
from typing import Any, Dict, List, Tuple, Union
import pyslha
import json
import itertools
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter

# Fix seed so results are reproducible!
np.random.seed(seed=123)



class effMap:
  def __init__(self, filepath: str) -> None:
    
    self.data = self.load_efficiency_map(filepath)
    self.eff_label = self.data.dtype.names[-1]
    self.vars_limits = self.setLimits(self.data)
    self.interp_nearest = self.setInterp(self.data)
    self.interp_hist2d = self.hist2d_lookup(self.data)
    self.interp_smooth = self.smooth_interpolator_2d(self.data, 
                                              sigma=0.0, 
                                              method='linear')
    
    
  def load_efficiency_map(self, csv_path: str) -> Any:
    """
    Load (pT, d0) -> efficiency points from the HEPData CSV and return
    an evaluator function: efficiency(pt, d0, clip=True).
    """

    with open(csv_path, 'r') as f:
      lines = [l for l in f.readlines() if not l.startswith('#')]
    data = np.genfromtxt(lines, delimiter=",", comments="#",
                        names=True,   dtype=float)
    return data
  
  @classmethod
  def centers_to_edges(cls,centers):
    centers = np.asarray(sorted(np.unique(centers)), dtype=float)
    if centers.size < 2:
        raise ValueError("Need at least two distinct centers to infer bin edges.")

    mids = 0.5 * (centers[:-1] + centers[1:])
    edges = np.empty(centers.size + 1, dtype=float)
    edges[1:-1] = mids
    edges[0] = centers[0] - 0.5 * (centers[1] - centers[0])
    edges[-1] = centers[-1] + 0.5 * (centers[-1] - centers[-2])

    return edges
  
  def setLimits(self,data):
    """
    Defines the allowed limits for computing efficiencies
    """

    if len(data.dtype.names) != 3:
      raise ValueError("data must have shape (N, 3): [x_center, y_center, content]")
    
    vars_limits = {}
    for var in data.dtype.names[:-1]:
      edges = effMap.centers_to_edges(data[var])
      vars_limits[var] = (edges[0], edges[-1])
  
    return vars_limits
  
  def hist2d_lookup(self, data):
    
    nbinsDict = {}
    edgesDict = {}
    for var in self.vars_limits.keys():
      nbinsDict[var] = len(data[var])
      edgesDict[var] = effMap.centers_to_edges(data[var])
    
    self.contents = np.full(tuple(list(nbinsDict.values())), 
                            np.nan, dtype=float)
    for row in data:
      x = row[data.dtype.names[0]]
      y = row[data.dtype.names[1]]
      val = row[self.eff_label]
      i = np.searchsorted(edgesDict[data.dtype.names[0]], x) - 1
      j = np.searchsorted(edgesDict[data.dtype.names[1]], y) - 1
      if 0 <= i < nbinsDict[data.dtype.names[0]] and 0 <= j < nbinsDict[data.dtype.names[1]]:
        self.contents[i, j] = val

    def get_bin_content(x, y, outside_value=np.nan):
        i = np.searchsorted(edgesDict[data.dtype.names[0]], 
                            x, side="right") - 1
        j = np.searchsorted(edgesDict[data.dtype.names[1]], y, side="right") - 1

        if i < 0 or i >= nbinsDict[data.dtype.names[0]] or j < 0 or j >= nbinsDict[data.dtype.names[1]]:
            return outside_value

        return self.contents[i, j]
    
    return get_bin_content
    
  
  def setInterp(self, data: Any) -> Any:

    interp = NearestNDInterpolator(list(zip(data[list(self.vars_limits.keys())[0]],
                                                 data[list(self.vars_limits.keys())[1]])),
                                                 data[self.eff_label])
    return interp

  def smooth_interpolator_2d(self, data, sigma=0.8, method='linear'):
    """
    Build a smooth interpolating function f(x, y) from 2D histogram bin contents.

    Parameters
    ----------
    xedges, yedges : array-like
        Bin edges for x and y with lengths Nx+1 and Ny+1.
    bin_values : array-like, shape (Nx, Ny)
        Bin contents defined on each 2D bin.
    sigma : float, optional
        Gaussian smoothing strength in bin units. Use 0 for no smoothing.
    method : {'linear', 'nearest'}, optional
        Interpolation mode for RegularGridInterpolator.
    fill_value : float or None, optional
        Value outside range. If None, performs extrapolation.

    Returns
    -------
    f : callable
        Function f(x, y) that accepts scalars or numpy arrays.
    """

    if len(data.dtype.names) != 3:
      raise ValueError("data must have shape (N, 3): [x_center, y_center, content]")
    
    # Define bin centers
    xcenters = np.unique(data[data.dtype.names[0]]).tolist()
    ycenters = np.unique(data[data.dtype.names[1]]).tolist()
    xedges = effMap.centers_to_edges(xcenters)
    yedges = effMap.centers_to_edges(ycenters)
    # The highest y edge (d0 = 300 mm) should have zero eff,
    # so add a point at last y edge so the interpolator interpolates to zero
    ycenters.append(yedges[-1])
    # The efficiency at very low y values (d0 < 25 mmm) grows much
    # faster than the linear behavior. Therefore we need to add a point
    # at the lowest y value with a larger efficiency to mimic this behavior
    ycenters.insert(0,yedges[0])

    mean_eff = np.zeros((len(xcenters), len(ycenters)))

    # Build a fast lookup from (pT, d0) -> row index in electron_reco.data
    centers = list(zip(data[data.dtype.names[0]], data[data.dtype.names[1]]))
    center_to_index = {center: idx for idx, center in enumerate(centers)}

    for i, j in itertools.product(range(len(xcenters)), range(len(ycenters))):
        idx = center_to_index.get((xcenters[i], ycenters[j]),None)
        if idx is None:
            continue
        mean_eff[i, j] = data[self.eff_label][idx]

    # Set the efficiency at the highest y value to zero
    # and a the lowest y value to 38% larger than the previous bin
    # to mimic the increase at low d0 
    # (the factor was estimated comparing the efficiencies from Fig.1 of 2011.07812 
    # to the efficiencies from FigAux_19a at the lowest d0 value)
    for i in range(mean_eff.shape[0]):
      mean_eff[i,-1] = 0.0
      mean_eff[i,0] = 1.5*mean_eff[i,1]

    
    z = mean_eff

    if z.shape != (len(xcenters), len(ycenters)):
        raise ValueError(
            f"bin_values shape {z.shape} is incompatible with edges "
            f"({len(xcenters)}, {len(ycenters)})."
        )


    # Fill NaNs before smoothing/interpolating to avoid artifacts.
    z_filled = np.nan_to_num(z, nan=0.0)

    if sigma and sigma > 0:
        z_smooth = gaussian_filter(z_filled, sigma=sigma, mode='nearest')
    else:
        z_smooth = z_filled


    interp = RegularGridInterpolator(
        (xcenters, ycenters),
        z_smooth,
        method=method,
        bounds_error=False,
        fill_value=None,
    )

    def f(x, y):
        x_arr, y_arr = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float))

        # Treat the histogram edge box as in-domain: clamp to nearest center
        # so points between edge and first/last center are not forced to fill_value.
        inside_edges = (
            (x_arr >= xedges[0]) & (x_arr <= xedges[-1]) &
            (y_arr >= yedges[0]) & (y_arr <= yedges[-1])
        )

        pts = np.column_stack([x_arr.ravel(), y_arr.ravel()])
        out = interp(pts).reshape(x_arr.shape)

        # Preserve fill_value behavior for points truly outside histogram edges.
        out = np.where(inside_edges, out, 0.0)

        return float(out) if out.ndim == 0 else out

    return f
  
  def efficiency(self,method='smooth',**kwargs) -> float:

    var_values = [kwargs.get(var,None) for var in self.vars_limits.keys()]
    if any(v is None for v in var_values):
      raise ValueError(f"Missing variable(s) for efficiency map: {self.vars_limits.keys()}. Got {kwargs.keys()}")
    if any(v < self.vars_limits[var][0] or v > self.vars_limits[var][1] for var,v in zip(self.vars_limits.keys(),var_values)):
      return 0.0
    if method == 'binned':
      eff = float(self.interp_hist2d(*var_values))
    elif method == 'nearest':
      eff = float(self.interp_nearest(*var_values))
    elif method == 'smooth':
      eff = float(self.interp_smooth(*var_values))
    else:
      raise ValueError(f"Method {method} not found!")
    if np.isnan(eff):
      return 0.0

    return eff



class cutFlow(object):

  def __init__(self,name: str,levels: List[str],zero_weight : Union[float,ndarray] = 0.0) -> None:
    self.name = name
    self.keys = levels[:]
    self.weights = np.array([zero_weight for _ in levels])
    self.weightsErr = np.array([zero_weight for _ in levels])
    self._current_level = 0

  def __repr__(self) -> str:
    return str(self)
  
  def __str__(self) -> str:
    return self.name

  def reset(self,to_level: int=0):
    self._current_level = to_level

  def fill_next(self,weight: Union[float,ndarray]):
    """
    Add weight to the next level in the cutflow
    """
    self._current_level += 1
    self.fill(weight)

  def fill_level(self,level: str,weight: Union[float,ndarray]):
    """
    Fill a specific level with the weight
    """
    if not (level in self.keys):
      raise ValueError(f"Level {level} not defined for {self}")
    clevel = self.keys.index(level)
    self.weights[clevel] += weight
    self.weightsErr[clevel] = np.sqrt(self.weightsErr[clevel]**2 + weight**2)
    
  def fill(self,weight: Union[float,ndarray]):
    """
    Add weight to the current level in the cutflow
    """
    
    clevel = self.current_level
    if not (0 <= clevel < len(self.weights)):
      msg = f"Trying to fill the cutflow at level {clevel+1}, but it only has {len(self.weights)} levels."
      msg += " Check your cutflow definitions and if it is being properly reset."
      raise ValueError(msg)
    self.weights[clevel] += weight
    self.weightsErr[clevel] = np.sqrt(self.weightsErr[clevel]**2 + weight**2)

  def divide(self,factor: float):
    """
    Divide cutflow levels by factor.
    """
    for iw,w in enumerate(self.weights):
      self.weights[iw] = w/factor
      self.weightsErr[iw] = self.weightsErr[iw]/factor

  @property
  def current_level(self) -> int:
    """
    Simple method for getting the current level of the cutflow
    """

    return self._current_level

  def to_dict(self) -> Dict[str, Union[Tuple[ndarray, ndarray], Tuple[float, float]]]:

    cDict = {k : (w,wErr) for k,w,wErr in zip(self.keys,self.weights,self.weightsErr)}

    return cDict
  
  def to_string(self) -> str:

    d = self.to_dict()
    lines = [f"==== {self.name} ==="]
    key_width = max(len(str(k)) for k in d) if d else 0
    for k,(w,wErr) in d.items():
      key_label = f"{k:<{key_width}}"
      if isinstance(w,(float,int)):
        lines.append(f"{key_label} = {w:1.4e} +- {wErr:1.4e}")
      elif isinstance(w,(list,ndarray)):
        l = f"{key_label} = "
        l += ' / '.join([f"{wx:1.4e} +- {wErrx:1.4e}" for wx,wErrx in zip(w,wErr)])
        lines.append(l)
    lines.append(f"===" + "="*len(self.name) + "===")
    return '\n'.join(lines)



#Object readers
def filterObjects(particleList: Any, 
                  pTmin: float, etaMax: float) -> List[Any]:
  filteredParticles = []
  for ptc in particleList:
    if ptc.PT < pTmin:
      continue
    if abs(ptc.Eta) > etaMax:
      continue

    filteredParticles.append(ptc)
  
  return filteredParticles


#Get Kinematic variables from objects
def deltaR(ptc1,ptc2) -> float:
  lv1 = ptc1.P4() #Check if this is the proper way to read 4vector and use as input for DeltaR
  lv2 = ptc2.P4()
  return lv1.DeltaR(lv2)

def DeltaPhi(ptc1, ptc2) -> float:
  return abs(ptc1.P4().DeltaPhi(ptc2.P4()))

def getD0(ptc) -> float:
  x = ptc.X
  y = ptc.Y
  phi = ptc.Phi
  pT = ptc.PT
  vTrack = np.array([x,y])
  pTrack = np.array([pT*np.cos(phi),pT*np.sin(phi)])
  d0 = np.linalg.norm(np.cross(vTrack,pTrack))/pT  
  return d0

def getZ0(ptc) -> float:
  x = ptc.X
  y = ptc.Y
  z = ptc.Z
  phi = ptc.Phi
  pT = ptc.PT
  pz = pT*np.sinh(ptc.Eta)
  vTrack = np.array([x,y])
  pTrack = np.array([pT*np.cos(phi),pT*np.sin(phi)])
  z0 = z - (pz/pT)*np.dot(vTrack,pTrack)/pT
  return z0

def getR(ptc) -> float:
  x = ptc.X
  y = ptc.Y
  return np.sqrt(x**2 + y**2)

def minDphilist(ptc1, listptc2, length, cut) -> float:
  if len(listptc2)==0:
    return 0
  infDphi = 99999999
  for iptc,ptc2 in enumerate(listptc2):
    if iptc>=length:
      break
    if ptc2.PT<cut:
      continue
    infDphi=min(infDphi,DeltaPhi(ptc1,ptc2))
  return infDphi

def overlapRemoval(input: List[Any],filter: List[Any],dR: float=0.4) -> List[Any]:
  
  if len(input)==0 or len(filter)==0:
    return input[:]
  
  output=[]
  for ptc1 in input:
    if any(deltaR(ptc1,ptc2)<dR for ptc2 in filter):
      continue
    output.append(ptc1)

  return output

def getLLPDecayRadius(llp) -> float:
  return np.sqrt(llp.daughter.X**2 + llp.daughter.Y**2)

def getLLPDecayTime(llp) -> float:
  return 1e9*(llp.daughter.T - llp.T) #Assume genpart T is in sec, convert to ns



def _iter_tracker_segments(layers, tol=1e-9):
    """
    Iterate over the layers and returns the layer boundaries as (rmin, rmax, zmin, zmax) tuples. For layers that are symmetric with respect to the z=0 plane, two segments are returned: one for the positive z side and one for the negative z side.
    """
    for layer in layers:
        rmin = layer['Rmin_mm']
        rmax = layer['Rmax_mm']
        zmin = layer['Zmin_mm']
        zmax = layer['Zmax_mm']

        if np.isclose(rmin, rmax, atol=tol) and np.isclose(zmin, 0.0, atol=tol):
            yield (rmin, rmax, -zmax, zmax)
            continue

        yield (rmin, rmax, zmin, zmax)
        if not np.isclose(zmin, 0.0, atol=tol) or not np.isclose(zmax, 0.0, atol=tol):
            yield (rmin, rmax, -zmax, -zmin)


def cross_layer(R0, z0, vR, vz, segment, tol=1e-9) -> bool:
    """
    Returns True if a track created at (R0, z0) and with velocity (vR, vz) crosses the layer defined by (rmin, rmax, zmin, zmax), and False otherwise.
    """
    rmin, rmax, zmin, zmax = segment

    # Check if is is an horizontal layer    
    if np.isclose(rmin, rmax, atol=tol):
        # A non-zero radial velocity is needed to cross it
        if np.isclose(vR, 0.0, atol=tol):
            return False

        # Compute the time to cross the layer's radius
        # (if it is negative, it means the track does not cross the layer)
        t_cross = (rmin - R0) / vR
        if t_cross <= tol:
            return False

        # Now compute the z coordinate at which the track crosses the layer, and check if it is within the layer's z range
        z_cross = z0 + vz * t_cross
        if (z_cross > zmax + tol) or (z_cross < zmin - tol):
            return False
        else:
            return True

    # Check if is is an horizontal layer   
    elif np.isclose(zmin, zmax, atol=tol):
        # A non-zero longitudinal velocity is needed to cross it
        if np.isclose(vz, 0.0, atol=tol):
            return False

        # Compute the time to cross the layer's radius
        # (if it is negative, it means the track does not cross the layer)
        t_cross = (zmin - z0) / vz
        if t_cross <= tol:
            return False

        # Now compute the radial coordinate at which the track crosses the layer, and check if it is within the layer's R range
        r_cross = R0 + vR * t_cross
        if (r_cross < rmin - tol) or (r_cross > rmax + tol):
            return False
        else:
            return True
    else:
        raise ValueError('Only axis-aligned tracker layers are supported.')


def count_crossed_layers(R0, z0, vR, vz, layers, tol=1e-9):
    """
    Returns the number of layers crossed by a track created at (R0, z0) and with velocity (vR, vz) among the given layers. The layers should be given as a list of (rmin, rmax, zmin, zmax) tuples, or as a structured array with fields 'Rmin_mm', 'Rmax_mm', 'Zmin_mm', 'Zmax_mm'. Only axis-aligned layers are supported.
    """
    ncrossings = 0
    for segment in _iter_tracker_segments(layers, tol=tol):
        cross = cross_layer(R0, z0, vR, vz, segment, tol=tol)
        if cross:
            ncrossings += 1
    return ncrossings


# ATLAS Pixel layers: 
pixelLayers = np.array([(33.5,33.5,0.,400.), (50.5,50.5,0.,400.), (88.5,88.5,0.,400.), (122.5,122.5,0.,400.),(88.5,130.,495.,495.),(88.5,130.,580.,580.),(88.5,130.,650.,650.)],dtype=np.dtype([('Rmin_mm', float),('Rmax_mm', float),('Zmin_mm', float),('Zmax_mm',float)]))
# ATLAS SCT layers: 
sctLayers = np.array([(299.0,299.0,0.,749.),(371.0,371.0,0.,749.),(443.0,443.0,0.,749.),(514.0,514.0,0.,749.),(337.6,560.0,853.8,853.8),(275.0,560.0,934.,934.0),(275.0,560.0,1091.5,1091.5),(275.0,560.0,1299.9,1299.9),(275.0,560.0,1399.7,1399.7),(275.0,560.0,1771.4,1771.4),(337.6,560.0,2115.2,2115.2),(408.0,560.0,2505.0,2505.0),(438.8,560.0,2720.2,2720.2)],dtype=np.dtype([('Rmin_mm', float),('Rmax_mm', float),('Zmin_mm', float),('Zmax_mm',float)]))

def count_tracker_layer_crossings(R0, z0, vR, vz, tol=1e-9):
    """
    Returns the number of pixel and SCT layers crossed by a track created at (R0, z0) and with velocity (vR, vz). The layers should be given as a list of (rmin, rmax, zmin, zmax) tuples, or as a structured array with fields 'Rmin_mm', 'Rmax_mm', 'Zmin_mm', 'Zmax_mm'. Only axis-aligned layers are supported.
    """

    pixel_hits = count_crossed_layers(R0, z0, vR, vz, pixelLayers, tol=tol)
    sct_hits = count_crossed_layers(R0, z0, vR, vz, sctLayers, tol=tol)
    return {
        'pixel': pixel_hits,
        'sct': sct_hits,
        'total': pixel_hits + sct_hits,
    }


def getModelInfo(bannerFile : str, llpPDG : int) -> Dict[str, Union[float,int]]:

    modelInfoDict : Dict[str, Union[float,int]] = {'llpPDG' : llpPDG}
    slhaData = None
    with open(bannerFile,'r') as ff:
        data = ff.read()
        if '<slha>' in data:
            slhaData = data.split('<slha>')[1].split('</slha>')[0]
            slhaData = pyslha.readSLHA(slhaData)
        if '<MGGenerationInfo>' in data:
            mgInfo = data.split('<MGGenerationInfo>')[1].split('</MGGenerationInfo>')[0]
            for l in mgInfo.split('\n'):
                l = l.strip()
                if not  l: continue
                k,v = l.split(':')
                modelInfoDict[k.replace('#','').strip()] = float(v)
        if '<MGRunCard>' in data:
            runInfo = data.split('<MGRunCard>')[1].split('</MGRunCard>')[0]
            fields = ['custom_fcts','pt_bias_target',
                      'pt_bias_enhancement_power', 'pt_bias_min']
            for l in runInfo.split('\n'):
                l = l.strip()
                if not  l: continue
                if l[0] == '#': #skip comments
                  continue
                for field in fields:
                  if field in l:
                    value = l.split('=')[0].strip()
                    if not value:
                      value = None
                    if field != 'custom_fcts' and value is not None:
                       value = float(value)
                    modelInfoDict[field] = value
                  
    if slhaData is None:
       raise ValueError(f'Error reading banner file {bannerFile}')
    
    modelInfoDict['mLLP'] = slhaData.blocks['MASS'][llpPDG]
    if slhaData.decays:
        modelInfoDict['tau0_ns'] = (6.58212e-16/slhaData.decays[llpPDG].totalwidth)
    
    return modelInfoDict

def saveOutput(resultsDict: Dict[str, Any],outputFile: str):

  saveDict = {}
  for k,v in resultsDict.items():
    if isinstance(v,np.ndarray):
      saveDict[k] = v.tolist()
    else:
       saveDict[k] = v
  # Sort results by type of key:
  saveDict = {k : v for k,v in sorted(list(saveDict.items()), key = lambda x: isinstance(x[1],list),
                                           reverse=False)}

        
  with open(outputFile, 'w') as f:
    json.dump(saveDict,f,indent=4)

def get_poi_upper_limit(signal_yields: Dict[str, float], 
                         bg_json_file : Path = Path('./ATLAS_data/likelihood/Comb_bkgonly.json'),
                         expected=None) -> float:

    import spey
    from spey_pyhf.helper_functions import WorkspaceInterpreter

    if expected is None:
       expected = spey.ExpectationType.observed

    pdf_wrapper = spey.get_backend('pyhf')
    with open(bg_json_file, "r") as f:
        bkg_only = json.load(f)
    interpreter = WorkspaceInterpreter(bkg_only)
    for channel in interpreter.channels:
        ns = signal_yields.get(channel,0.)
        interpreter.inject_signal(channel, [ns])
    statistical_model = pdf_wrapper(
        analysis="simple_pyhf",
        background_only_model=interpreter.background_only_model,
        signal_patch=interpreter.make_patch(),
    )
    mu_UL = statistical_model.poi_upper_limit(expected=expected)

    return mu_UL