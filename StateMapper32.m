classdef StateMapper32 < handle
% This class facilitates working with graphs that contain states representable by 32-bit numbers.

properties (Access = private, Constant = true)
  NBITS = int8(32);
  UINT_CLASS = 'uint32';
  UINT_CLASS_FUNC = @uint32;
end

properties (GetAccess = public, SetAccess = immutable)
  configBits   (1,1) int8
  timestepBits (1,1) int8
end

properties (GetAccess = private, SetAccess = immutable)
  bitmaskStateFrom (1,1) uint32
  bitmaskStateTo   (1,1) uint32
  bitmaskTimestep  (1,1) uint32 = intmax('uint32');  
end

methods
  function sm = StateMapper32(configBits, timestepBits)
    arguments
      configBits   (1,1) int8 = 8;  % up-to 256 configurations
      timestepBits (1,1) int8 = 13; % up-to 8192 steps
    end
    % Input testing:
    spareBits = sm.NBITS - (2*configBits + timestepBits);
    assert( spareBits >= 0, 'StateMapper32:tooManyBits', ...
      'The specified number of states cannot fit into %d bits!', sm.NBITS);
    % Store settings:
    sm.configBits = configBits;
    sm.timestepBits = timestepBits;
    % Compute masks:
    sm.bitmaskStateTo   = 2^sm.UINT_CLASS_FUNC(sm.NBITS-spareBits-timestepBits)           -1; % dec2bin(sm.bitmaskStateTo,  sm.NBITS)
    sm.bitmaskStateFrom = 2^sm.UINT_CLASS_FUNC(sm.NBITS-spareBits-timestepBits-configBits)-1; % dec2bin(sm.bitmaskStateFrom,sm.NBITS)
  end
end

methods (Access = public)
  function uintState = toUint(sm, stateFrom, stateTo, timeStamp)
    %% Input tests:
    assert(all(stateFrom <= p2(sm.configBits)) ...
        && all( stateTo  <= p2(sm.configBits)),...
      'toUint:inputStateOverflow', ...
      'One or more states IDs cannot be represented by the available number of bits.');
    assert(all(timeStamp <= p2(sm.timestepBits)),'toUint:inputTimestepOverflow',...
      'One or more timestamps cannot be represented by the available number of bits.');
    assert(isequal(numel(stateFrom), numel(stateTo), numel(timeStamp)), ...
      'toUint:differentInputLengths',...
      'All inputs must contain the same number of elements!');
    %% Actual computation:
    uintState = toUintUnsafe(sm, stateFrom, stateTo, timeStamp);
  end
  
  function uintState = toUintUnsafe(sm, stateFrom, stateTo, timeStamp)
    % The business logic of toUint(). Can be invoked directly (should be slightly faster)  
    % but has the risk of producing nonsense.      
    uintState = ...
                sm.UINT_CLASS_FUNC(stateFrom)                   + ...
      bitshift( sm.UINT_CLASS_FUNC( stateTo ),   sm.configBits) + ...
      bitshift( sm.UINT_CLASS_FUNC(timeStamp), 2*sm.configBits);  
  end
  
  function [stateFrom, stateTo, timeStamp] = fromUint(sm, uintState)
    arguments
      sm (1,1) StateMapper32
      uintState (:,:) uint32 % dec2bin(uintState,32)
    end
    timeStamp = bitshift( bitand( sm.bitmaskTimestep,  uintState, sm.UINT_CLASS), -2*sm.configBits);
    stateTo   = bitshift( bitand( sm.bitmaskStateTo,   uintState, sm.UINT_CLASS), -1*sm.configBits);
    stateFrom = bitshift( bitand( sm.bitmaskStateFrom, uintState, sm.UINT_CLASS), -0*sm.configBits);
  end
end

end

function out = p2(in)
  out = bitshift(1, in, StateMapper32.UINT_CLASS);
end
