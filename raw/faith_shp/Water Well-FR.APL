<?xml version="1.0" encoding="UTF-8"?>
<ArcPad>
   <LAYER name="water well">
      <SYMBOLOGY>
         <SIMPLELABELRENDERER visible="false" field="WELL_NAME">
            <TEXTSYMBOL angle="0" fontcolor="0,0,0" font="Arial" fontsize="8" horzalignment="left" vertalignment="bottom" rtl="false">
            </TEXTSYMBOL>
         </SIMPLELABELRENDERER>
         <VALUEMAPRENDERER lookupfield="Status">
            <EXACT value="Abandoned" label="Abandoned">
               <GROUPSYMBOL>
                  <TRUETYPEMARKERSYMBOL font="ESRI Geometric Symbols" character="89" fontcolor="78,78,78" fontsize="13" angle="0">
                  </TRUETYPEMARKERSYMBOL>
                  <TRUETYPEMARKERSYMBOL font="ESRI Geometric Symbols" character="207" fontcolor="255,0,197" fontsize="13" angle="0">
                  </TRUETYPEMARKERSYMBOL>
               </GROUPSYMBOL>
            </EXACT>
            <EXACT value="Non Producing" label="Non Producing">
               <GROUPSYMBOL>
                  <TRUETYPEMARKERSYMBOL font="ESRI Geometric Symbols" character="89" fontcolor="78,78,78" fontsize="14" angle="0">
                  </TRUETYPEMARKERSYMBOL>
                  <TRUETYPEMARKERSYMBOL font="ESRI Geometric Symbols" character="35" fontcolor="252,196,179" fontsize="14" angle="0">
                  </TRUETYPEMARKERSYMBOL>
               </GROUPSYMBOL>
            </EXACT>
            <EXACT value="Producing" label="Producing">
               <GROUPSYMBOL>
                  <TRUETYPEMARKERSYMBOL font="ESRI Geometric Symbols" character="89" fontcolor="78,78,78" fontsize="18" angle="0">
                  </TRUETYPEMARKERSYMBOL>
                  <TRUETYPEMARKERSYMBOL font="ESRI Geometric Symbols" character="35" fontcolor="197,0,255" fontsize="18" angle="0">
                  </TRUETYPEMARKERSYMBOL>
               </GROUPSYMBOL>
            </EXACT>
            <OTHER label="Working Water Well">
               <GROUPSYMBOL>
                  <TRUETYPEMARKERSYMBOL font="ESRI Geometric Symbols" character="89" fontcolor="78,78,78" fontsize="17" angle="0">
                  </TRUETYPEMARKERSYMBOL>
                  <TRUETYPEMARKERSYMBOL font="ESRI Geometric Symbols" character="35" fontcolor="196,149,254" fontsize="17" angle="0">
                  </TRUETYPEMARKERSYMBOL>
               </GROUPSYMBOL>
            </OTHER>
         </VALUEMAPRENDERER>
      </SYMBOLOGY>
   </LAYER>
</ArcPad>
