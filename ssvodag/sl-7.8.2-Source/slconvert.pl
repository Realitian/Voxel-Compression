#!/bin/perl
# +++HDR+++
# ======================================================================
#   This file is part of the SL software library.
# 
#     Copyright (C) 1993-2002 by Enrico Gobbetti (gobbetti@crs4.it)
#     Copyright (C) 1996-2002 by CRS4, Cagliari, Italy.
# 
#   For more information, visit the CRS4 Visualization and Virtual 
#   Reality Group web pages at http://www.crs4.it/vvr/.
# ======================================================================
# ---HDR---

$sedscript = "" 
. " -e 's#\"slTuple\.h\"#<sl/fixed_size_array.hpp>#g'" 
. " -e 's#slTuple2i#sl::fixed_size_array<2,int>#g'"
. " -e 's#slTuple3i#sl::fixed_size_array<3,int>#g'"
. " -e 's#slTuple4i#sl::fixed_size_array<4,int>#g'"
. " -e 's#sl::fixed_size_array<2,int>#sl::tuple2i#g'"
. " -e 's#sl::fixed_size_array<3,int>#sl::tuple3i#g'"
. " -e 's#sl::fixed_size_array<4,int>#sl::tuple4i#g'"
. " -e 's#\"slPoint\.h\"#<sl/fixed_size_point.hpp>#g'" 
. " -e 's#slPoint2f#sl::point2f#g'" 
. " -e 's#slPoint3f#sl::point3f#g'" 
. " -e 's#slColor3f#sl::color3f#g'" 
. " -e 's#slColor4f#sl::color4f#g'" 
. " -e 's#\"slVector\.h\"#<sl/fixed_size_vector.hpp>#g'" 
. " -e 's#slVector2f#sl::column_vector2f#g'" 
. " -e 's#slVector3f#sl::column_vector3f#g'" 
. " -e 's#\"slFrame\.h\"#<sl/affine_map.hpp>#g'" 
. " -e 's#slFrame3f#sl::affine_map3f#g'" 
. " -e 's#sl::affine_map3f::#sl::affine_map_factory3f::#g'" 
. " -e 's#\"slMatrix\.h\"#<sl/projective_map.hpp>#g'" 
. " -e 's#slMatrix4f#sl::projective_map3f#g'" 
. " -e 's#\"slBox\.h\"#<sl/axis_aligned_box.hpp>#g'" 
. " -e 's#\slBox3f#sl::aabox3f#g'"
. " -e 's#box\.size#box.diagonal#g'" 
. " -e 's#toZero#to_zero#g'" 
. " -e 's#toEmpty#to_empty#g'" 
. " -e 's#toPointer#to_pointer#g'" 
. " -e 's#toIdentity#to_identity#g'" 
. " -e 's#twoNorm#two_norm#g'" 
. " -e 's#oneNorm#one_norm#g'" 
. " -e 's#transformedBy#transformed_by#g'"
. " -e 's#\<iostream.h\>#\<iostream\>#g'" 
. " -e 's#cerr#std::cerr#g'" 
. " -e 's#cout#std::cout#g'" 
. " -e 's#endl#std::endl#g'" 
. " -e 's#\<vector.h\>#\<vector\>#g'"
. " -e 's#vector\<#std::vector\<#g'"
. " -e 's#ostream#std::ostream#g'"
. " -e 's#istd::ostream#iostream#g'"
. " -e 's#pair\<#std::pair\<#g'"
. " -e 's#make_pair#std::make_pair#g'"
. " -e 's#std::size_t#size_t#g'"
. " -e 's#std::ptrdiff_t#ptrdiff_t#g'"
. " -e 's#::std::#::#g'"
. " -e 's#_std::#_#g'"
. " -e 's#std::std::#std::#g'"
;


sub verbose {
  print STDERR "* @_\n";
}

sub sysExec {
#    verbose "Command = @_\n";
    system(@_);
}

sub process {
    while (<*.{h,hpp,cpp,cc}>) {
        my $fin =  $_;
        my $fbak = $fin.".bak";
	my $fout = $fin;

        verbose "Backup: $fin --> $fbak";
        sysExec("cp $fin $fbak");
	verbose("converting: $fin -> $fout");
	sysExec("cat $fbak | sed $sedscript > $fout"); 
    }
}

process;
  
