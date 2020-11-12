classdef Season < uint8
  enumeration
    Winter  (1)
    Spring  (2)
    Summer  (3)
    Autumn  (4)
  end
  
  methods (Access = public, Static = true)
    function s = whichSeason(timestamp, hemisphere)
      arguments
        timestamp (:,:) datetime
        hemisphere (1,1) char {mustBeMember(hemisphere, {'N','S'})} = 'N'
      end
      N_SEASONS = 4;
      N_MONTHS = 12;
      N_MONTHS_PER_SEASON = N_MONTHS / N_SEASONS;
      sz = size(timestamp);
      mo = month(timestamp(:)) == 1:N_MONTHS;
      [~,I] = max(any(reshape(circshift(mo, 1, 2), [], N_MONTHS_PER_SEASON, N_SEASONS),2),[],3);
      % Explanation of the above line:
      % 1) circshift is applied so that all winter months are together [1, 2 ,... 12] -> [12, 1, 2, ...]
      % 2) the result is reshaped into a 3d array where each slice represents a season
      % 3) for each timestamp, it is checked whether any month of a season contains a true value
      % 4) max is used to find the slice in which a true value was found
      switch hemisphere
        case 'N' 
          s = reshape( Season(I), sz);
        case 'S'
          s = reshape( Season( mod(I-3,4)+1 ), sz);
      end
    end
    
    function tf = isSummer(timestamp, hemisphere)
      arguments
        timestamp (:,:) datetime
        hemisphere (1,1) char {mustBeMember(hemisphere, {'N','S'})} = 'N'
      end
      sz = size(timestamp);
      switch hemisphere
        case 'N'          
          tf = Season.isSummerNH(timestamp(:));
        case 'S'
          tf = Season.isSummerSH(timestamp(:));
      end
      tf = reshape(tf, sz);
    end
    
    function tf = isAutumn(timestamp, hemisphere)
      arguments
        timestamp (:,:) datetime
        hemisphere (1,1) char {mustBeMember(hemisphere, {'N','S'})} = 'N'
      end
      sz = size(timestamp);
      switch hemisphere
        case 'N'          
          tf = Season.isAutumnNH(timestamp(:));
        case 'S'
          tf = Season.isAutumnSH(timestamp(:));
      end
      tf = reshape(tf, sz);
    end
    
    function tf = isWinter(timestamp, hemisphere)
      arguments
        timestamp (:,:) datetime
        hemisphere (1,1) char {mustBeMember(hemisphere, {'N','S'})} = 'N'
      end
      sz = size(timestamp);
      switch hemisphere
        case 'N'          
          tf = Season.isWinterNH(timestamp(:));
        case 'S'
          tf = Season.isWinterSH(timestamp(:));
      end
      tf = reshape(tf, sz);
    end
    
    function tf = isSpring(timestamp, hemisphere)
      arguments
        timestamp (:,:) datetime
        hemisphere (1,1) char {mustBeMember(hemisphere, {'N','S'})} = 'N'
      end
      sz = size(timestamp);
      switch hemisphere
        case 'N'          
          tf = Season.isSpringNH(timestamp(:));
        case 'S'
          tf = Season.isSpringSH(timestamp(:));
      end
      tf = reshape(tf, sz);
    end
    
  end
  
  methods (Access = private, Static = true)
    %% Northern hemisphere:
    function tf = isSummerNH(timestamp)
      mo = month(timestamp);
      tf = any(mo == 6:8);
    end
    function tf = isAutumnNH(timestamp)
      mo = month(timestamp);
      tf = any(mo == 9:11);
    end
    function tf = isWinterNH(timestamp)
      mo = month(timestamp);
      tf = any(mo == [12,1,2]);
    end
    function tf = isSpringNH(timestamp)
      mo = month(timestamp);
      tf = any(mo == 3:5);
    end
    %% Southern hemisphere:
    function tf = isSummerSH(timestamp)
      tf = Season.isWinterNH(timestamp);
    end
    function tf = isAutumnSH(timestamp)
      tf = Season.isSpringNH(timestamp);
    end
    function tf = isWinterSH(timestamp)
      tf = Season.isSummerNH(timestamp);
    end
    function tf = isSpringSH(timestamp)
      tf = Season.isAutumnNH(timestamp);
    end
  end
end
  