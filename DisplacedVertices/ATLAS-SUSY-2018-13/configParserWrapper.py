#!/usr/bin/env python3

#A wrapper for ConfigParser which allows to evaluate expressions
#in the parameters. The expressions should be of the type ${expr}.
#References for other parameters in the parser should be in the format section:parameter.


from math import *
import re, itertools
import numpy
import logging
logger = logging.getLogger("MG5Scan")

try:
    from ConfigParser import RawConfigParser,InterpolationDepthError,ParsingError
except:
    from configparser import RawConfigParser,InterpolationDepthError,ParsingError


class ConfigParserExt(RawConfigParser):
    
    
    def __init__(self,*args,**kargs):
        RawConfigParser.__init__(self,*args,**kargs)
        self.cur_depth = 0
        self.MAX_INTERPOLATION_DEPTH=100
        self.optionxform=str    #Preserve string cases    

    def toDict(self,raw=True):
        """
        Convert parser to dictionary.
        If raw=True, will return the raw values, otherwise it will
        try to evaluate them.
        """
        
        parserDict = {}
        for s in self.sections():
            parserDict[s] = {}
            for var in self.options(s):
                parserDict[s][var] = self.get(s,var,raw=raw)
        
        return parserDict
    
    def read_dict(self,parserDict):
        """
        Uses the dictionary to set the sections, options and option
        values. parserDict should be a dictionary in the format {'section1' : {'option1' : value1, 'option2' : value2,...}}
         
        :parameter parserDict: Dictionary with values for the sections and options.
         
        """
         
        for section in parserDict:
            if not section in self.sections():
                self.add_section(section)
            for option,value in parserDict[section].items():
                if isinstance(value,str):
                    self.set(section,option,value)
                else:
                    self.set(section,option,str(value))

    
    def get(self, section, option, raw=False, current_depth=0):
        
        valueRaw = RawConfigParser.get(self, section, option)
        if raw:            
            return valueRaw

        self.cur_depth = current_depth

        #Find all instances of ${(alphanumeric characters):(alphanumeric characters)}
        ret = valueRaw
        varSectOpts = re.findall(r'\$\{(\w*):(\w*)\}',ret)
        if varSectOpts:
            for v_section, v_option in varSectOpts:
                self.cur_depth = self.cur_depth + 1
                #Skip matches which do not correspond to any section
                if not v_section in self.sections():
                    continue
                #Skip matches which do not correspond to any option
                if not v_option in self.options(v_section):
                    continue
                
                #Try to get the value for the corresponding section/option
                #(avoid infinite loops with interpolation depth:            
                self.cur_depth = self.cur_depth + 1 
                if self.cur_depth < self.MAX_INTERPOLATION_DEPTH:
                    sub = self.get(v_section, v_option, current_depth=self.cur_depth)
                    #Replace all instances of (section):(option)
                    ret = re.sub(r'\$\{%s:%s\}' %(v_section, v_option), str(sub),ret)
                else:
                    raise InterpolationDepthError(option, section, valueRaw)
                
        #Now try to get all variables from the same section:
        for v_option in self.options(section):            
            #Check if the option appears in the string surrounded by any non-alphanumeric character
            #or the beginning of the expression or end of the expression
            if re.findall(r'\$\{%s\}'%v_option,ret):
                self.cur_depth = self.cur_depth + 1 
                if self.cur_depth < self.MAX_INTERPOLATION_DEPTH:
                    sub = self.get(section, v_option, current_depth=self.cur_depth)
                    #Replaces the correct instances:
                    ret = re.sub(r'\$\{%s\}' %v_option, str(sub),ret)
                else:
                    raise InterpolationDepthError(option, section, valueRaw)
        
        
        try:
            return eval(ret)
        except:
            return ret

        
    def getstr(self,*args,**kargs):
        
        return str(self.get(*args,**kargs))
    
    
    def expandLoops(self):
        """
        If one or more options have the tag $loop{}, it will
        generate new parsers for product of all values in the loop.
        E.g. if option1 = $loop{[100,200]} and option2 = $loop{['a','b']},
        it will generate parsers with option1,option2 = (100,'a'), (100,'b'), (200,'a') and (200,'b')
        
        :return: List of parsers with the options set to the respective value defined by the loop
        """
        
        loopVars = []
        varValues = []
        #Collect all loop options
        for section in self.sections():
            for option in self.options(section):
                ret = self.get(section,option,raw=True)
                varSectLoop = re.findall(r'\$loop\{(.*)\}',ret)
                if not varSectLoop:
                    continue
                if len(varSectLoop) > 1:
                    raise ParsingError("Syntax error. Multiple loops found for option %s." %option)
                loopStr = varSectLoop[0]
                try:
                    varList = eval(loopStr)
                except:
                    raise ParsingError("Could not evaluate loop %s" %loopStr)
                if not isinstance(varList,(list,numpy.ndarray,tuple)):
                    raise ParsingError("Loop expression %s did not generate a list" %loopStr)
                loopVars.append((section,option))
                varValues.append(varList)
                
        if not loopVars:
            return [self]

        logger.info(" Looping over variables:  " + ", ".join(["%s:%s" %lvar for lvar in loopVars]))
        parserList = []
        for values in itertools.product(*varValues):
            newParser = ConfigParserExt()
            newParser.read_dict(self.toDict(raw=True))
            for i,v in enumerate(values):
                sect,opt = loopVars[i]   
                newParser.set(sect,opt,str(v))
            parserList.append(newParser)
            
        return parserList

if __name__ == "__main__":
    
    parFile = './test/test.ini'
    parFile = './test.ini'
    
    parser = ConfigParserExt()
    ret = parser.read(parFile)
#     ret = parser.read('slha-parameters.ini')
    
    print(parser.get("MadGraphPars","mg5out"))
#     print(parser.sections())
    print(parser.get('slhaCreator','inputFile'))
#     print(parser.get('slhaCreator','slhaout'))
#     print(parser.get('MadGraphSet','F'))
#     print(parser.get('MadGraphSet','MH'))
#     print(parser.get('MadGraphSet','v'))
#     print(parser.get('MadGraphSet','MHu'))
#     print(parser.get('MadGraphSet','MHTodd'))
#     print(parser.get('MadGraphSet','MT'))
#     
    
            
