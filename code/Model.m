classdef Model < handle
    %MODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        M
        C
        K
        F
        input_file
        npos

    end
    properties (Access = private)
        
    end
    methods
        function obj = Model(input_filename)
            if nargin == 1
                obj.Assem(input_filename);
                obj.input_file = input_filename;
            end
        end
    end
    
end

