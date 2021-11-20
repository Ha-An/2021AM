################################################################################################################### 
# Creator: Yosep Oh (yosep.oh@kgu.ac.kr) 
# Environment: Blender 2.82; Windows 10 Enterprise; 
# Requirement (add-on script): V-HACD (https://github.com/kmammou/v-hacd)
# Description: This code is a python script for Blender 2.82. This script represents the evaluation indicators of part decomposition for Additive Manufacturing. 
################################################################################################################### 
import bpy, bmesh 
from math import radians, degrees, cos, inf
from scipy.spatial import distance
from statistics import mean
from mathutils import Vector 
import mathutils
import math
import timeit
import random as rd

# Setting parameters  
THRESHOLD_SURFACE_ROUGHNESS = 45    # degrees
THRESHOLD_SHARP_EDGE = 170          # degrees
THRESHOLD_OVERHANG = 60             # degrees 
THRESHOLD_GAP_DISTANCE = 2           
MAX_PART_SIZE = [200, 200, 300]     
THRESHOLD_MIN_SIZE = 1              
# Computation time
tic = timeit.default_timer() 

# ROUHGNESS #####################################################################################################
def roughness_for_each_part(part, CH, buildOrientationNormal):     
    # Calculate surface roughness
    bmPart = bmesh.new()
    bmPart.from_mesh(part.data) 
    surfaceRoughness =0
    bmPart.faces.ensure_lookup_table() 
    for f in bmPart.faces:   
        global_f_normal = f.normal
        if global_f_normal != Vector((0.0, 0.0, 0.0)):
            theta = buildOrientationNormal.angle(global_f_normal)    
            if 0<theta<radians(180):    
                SR = abs(cos(theta))
            else: 
                SR = 0 
            surfaceRoughness += SR * f.calc_area()    
    bmPart.free()   
    return surfaceRoughness  
def roughness_of_part_surface(decomposedParts): 
    sum=0    
    for i in decomposedParts:     
        sum += roughness_for_each_part(i["part"], i["CH"], i["buildOrientationNormal"])
    return sum 

# INTENSITY (OVERHANG)  ##########################################################################################
def intensity_of_critical_facets_for_each_part(part, CH, buildOrientationNormal):  
    # Calculate overhang area
    bmPart = bmesh.new()
    bmPart.from_mesh(part.data) 
    unsupportedSurfaceIntensity =0 
    bmPart.faces.ensure_lookup_table()   
    for f in bmPart.faces:   
        global_f_normal = f.normal
        if global_f_normal != Vector((0.0, 0.0, 0.0)):
            theta = buildOrientationNormal.angle(global_f_normal)   
            if radians(90) < theta:
                if radians(180 - THRESHOLD_OVERHANG) < theta and theta < radians(180-0.05):  ## 0.5 is an adjustment value for implementation
                    OH = abs(cos(theta))   
                else: 
                    OH = 0 
                unsupportedSurfaceIntensity += OH * f.calc_area()     
    bmPart.free()  
    return unsupportedSurfaceIntensity
# Sum of the area of causing overhangs for all convex hulls
def intensity_of_critical_facets(decomposedParts): 
    sum=0    
    for i in decomposedParts:     
        a = intensity_of_critical_facets_for_each_part(i["part"], i["CH"], i["buildOrientationNormal"])
        sum += a
    return sum 

# SHARPNESS ################################################################################################
def length_of_sharp_edges_for_each_part(part):      
    bmPart = bmesh.new()
    bmPart.from_mesh(part.data) 
    sharpEdgeLength =0
    bmPart.edges.ensure_lookup_table() 
    for e in bmPart.edges:  
        f1, f2 = e.link_faces
        if f1.normal != Vector((0,0,0)) and f2.normal != Vector((0,0,0)):
            theta = degrees(f1.normal.angle(f2.normal)) 
            if theta > THRESHOLD_SHARP_EDGE:
                sharpEdgeLength += e.calc_length()    
    bmPart.free()   
    return sharpEdgeLength  
def length_of_sharp_edges(decomposedParts): 
    sum=0    
    for i in decomposedParts:     
        sum += length_of_sharp_edges_for_each_part(i["part"])
    return sum 

# GAP  ##################################################################################################
def area_of_facets_causing_critical_gaps_for_each_part(part):  
    mat_global_obj = part.matrix_world  
    # Creating mesh data for the input part
    bmObj = bmesh.new()
    bmObj.from_mesh(part.data)
    bmObj.faces.ensure_lookup_table()   
    nonAllowableGapArea = 0 
    for f in bmObj.faces:  
        # distance check by using RayCast 
        c = f.calc_center_median()
        p = c +f.normal*0.001 
        hit_data = part.ray_cast(p, f.normal)      
        if hit_data[0]:    
            gapDistance = distance.euclidean(c, hit_data[1])  
            if gapDistance < THRESHOLD_GAP_DISTANCE: 
                nonAllowableGapArea +=f.calc_area()     
    bmObj.free()                
    return nonAllowableGapArea 
def area_of_facets_causing_critical_gaps(decomposedParts): 
    sum=0    
    for i in decomposedParts:    
        a = area_of_facets_causing_critical_gaps_for_each_part(i["part"])
        sum += a
    return sum  

# CONCAVITY #########################################################################################################
def volume_of_concave_spaces_for_each_part(part, CH): 
    bmPart = bmesh.new()
    bmPart.from_mesh(part.data)   
    volumePart = bmPart.calc_volume(signed=True)  
    bmCH = bmesh.new()
    bmCH.from_mesh(CH.data)   
    volumeCH = bmCH.calc_volume(signed=True)   
    return max((volumeCH - volumePart), 0) 
# Sum of concave space volume for all convex hulls
def volume_of_concave_spaces(decomposedParts): 
    sum=0    
    for i in decomposedParts:        
        v=volume_of_concave_spaces_for_each_part(i["part"], i["CH"])
        sum += v
    return sum 
 
# INTERFACE #################################################################################################
def area_of_connection_interface(decomposedParts): 
    sum=0    
    if len(decomposedParts) >1:
        for i in decomposedParts:         
            sum +=i["connectionArea"]
        contactArea = sum / 2 # Two parts share a contact area
    else:
        contactArea = 0
    return contactArea  
# FEASIBILITY ###############################################################################################
def feasibility_of_decomposed_parts(decomposedParts):  
    M = 0
    for i in decomposedParts:      
        maxPartDim = max(i["part"].dimensions) 
        minWorkspaceDim = min(MAX_PART_SIZE)       
        if maxPartDim < minWorkspaceDim:
            M = 0
        else:
            M = inf 
            print("i: ", i["part"])   
            print("maxPartDim: ", maxPartDim)      
            print("minWorkspaceDim: ", minWorkspaceDim)    
            break 
    return M  

# Part decomposition for convex features ######################################################################################
def bisectMaxP(maxP, dir, planeLocation, tmp_planeNormal): 
    if dir == "top":
        planeNormal = tmp_planeNormal
    else:
        planeNormal = -tmp_planeNormal 
    bpy.ops.object.select_all(action='DESELECT')  
    maxP.select_set(True) 
    bpy.ops.object.duplicate()     
    part=bpy.context.selected_objects[0]   
    bpy.context.view_layer.objects.active = part 
    bpy.ops.object.mode_set(mode='EDIT')  
    bpy.ops.mesh.select_all(action='SELECT')  
    bpy.ops.mesh.bisect(plane_co=planeLocation,plane_no=planeNormal,clear_outer=False, use_fill=True, clear_inner=True) 
    bpy.ops.object.mode_set(mode='OBJECT')
    # Getting the area of connection interface
    if dir == "top": 
        ray_begin = planeLocation - planeNormal * 0.001
        ray_begin_local = part.matrix_world.inverted() @ ray_begin  # covert ray_begin to "plane_ob" local space
        hit_data = part.ray_cast(ray_begin_local, planeNormal) 
        if hit_data[0]:
            connectionFaceIndex = hit_data[3]   
            bm = bmesh.new()
            bm.from_mesh(part.data)   
            bm.faces.ensure_lookup_table()
            for f in bm.faces: 
                if f.index == connectionFaceIndex:
                    connectionArea = f.calc_area()
            bm.free() 
        else:
            connectionArea=0 
    # Getting the convex hull of the part    
    bpy.ops.object.duplicate()     
    partCH=bpy.context.selected_objects[0]   
    partCH.name=partCH.name+"_hull"
    bpy.context.view_layer.objects.active = partCH 
    bpy.ops.object.mode_set(mode='EDIT')  
    bpy.ops.mesh.select_all(action='SELECT')
    bpy.ops.mesh.convex_hull()
    bpy.ops.object.mode_set(mode='OBJECT')
    bpy.context.object.display_type="WIRE"
    # Getting the buildOrientationNormal of the part    
    bmCH = bmesh.new()
    bmCH.from_mesh(partCH.data)   
    bmCH.faces.ensure_lookup_table()
    maxFace = 0
    buildOrientationNormal = Vector((0.0, 0.0, 0.0))
    for f in bmCH.faces:
        if f.calc_area() > maxFace:
            maxFace = f.calc_area()
            buildOrientationNormal = -f.normal.copy() 
    bmCH.free()  
    # Triangulating the decomposed parts
    bm = bmesh.new()
    bm.from_mesh(partCH.data) 
    bmesh.ops.triangulate(bm, faces=bm.faces[:], quad_method='BEAUTY', ngon_method='BEAUTY')
    bm.to_mesh(partCH.data)
    bm.free()  
    bm = bmesh.new()
    bm.from_mesh(part.data) 
    bmesh.ops.triangulate(bm, faces=bm.faces[:], quad_method='BEAUTY', ngon_method='BEAUTY')
    bm.to_mesh(part.data)
    bm.free() 
    # Changing object color 
    bpy.ops.object.select_all(action='DESELECT')  
    part.select_set(True)
    bpy.context.view_layer.objects.active = part   
    bpy.context.object.color = (rd.random(), rd.random(), rd.random(), 1) 
    bpy.ops.object.origin_set( type = 'ORIGIN_GEOMETRY' ) 
    if dir == "top":
        return part, partCH, buildOrientationNormal, connectionArea
    else: 
        return part, partCH, buildOrientationNormal
def plane_normal(maxPart):
    (Xd, Yd, Zd) = maxPart.dimensions
    XY = Xd*Yd
    YZ = Yd*Zd
    XZ = Xd*Zd               
    if min(XY, YZ, XZ) == YZ:
        planeNormal = Vector((-1, 0, 0)).copy()
    elif min(XY, YZ, XZ) == XZ:
        planeNormal = Vector((0, -1, 0)).copy()
    elif min(XY, YZ, XZ) == XY:
        planeNormal = Vector((0, 0, -1)).copy()
    return planeNormal 
def PD_for_convex_features(input_parts):
    decomposedParts=[]
    manifoldAll=True
     
    for i in input_parts:  
        if max(i["part"].dimensions) > min(MAX_PART_SIZE):  
            convexParts=[{"part": i["part"], "CH": i["CH"], "buildOrientationNormal": i["buildOrientationNormal"], "connectionArea": i["connectionArea"]}]
            largeParts=True 
            planeLocation = convexParts[0]["part"].location    
            manifoldCount =0
            
            maxP=convexParts[0]  
            planeNormal = plane_normal(maxP["part"])
             
            while(largeParts==True and manifoldCount <=20):
                print(" ")
                print("manifoldCount ", manifoldCount)   
                print("maxP ", maxP["part"].name) 
                print("dimensions ", maxP["part"].dimensions)   
                       
                # Getting the top part   
                print("planeLocation ", planeLocation) 
                print("planeNormal ", planeNormal) 
                topPart, topPartCH, topBuildOrientationNormal, connectionArea = bisectMaxP(maxP["part"], "top", planeLocation, planeNormal) 
                # Getting the bottom part   
                botPart, botPartCH, botBuildOrientationNormal = bisectMaxP(maxP["part"], "bottom", planeLocation, planeNormal)
                print("topPart ", topPart.dimensions) 
                print("botPart ", botPart.dimensions) 
                # Checking a manifold geometry
                bmTop = bmesh.new()
                bmTop.from_mesh(topPart.data)   
                bmTop.edges.ensure_lookup_table()
                bmBot = bmesh.new()
                bmBot.from_mesh(botPart.data)   
                bmBot.edges.ensure_lookup_table() 
                manifoldAll=True
                for e in bmTop.edges:
                    if e.is_manifold==False:
                        manifoldAll=False 
                for e in bmBot.edges:
                    if e.is_manifold==False:
                        manifoldAll=False 
                if manifoldAll == True:
                    connectionAreaSum = connectionArea + i["connectionArea"]
                    convexParts.append({"part":topPart, "CH":topPartCH, "buildOrientationNormal": topBuildOrientationNormal, "connectionArea": connectionAreaSum}) 
                    convexParts.append({"part":botPart, "CH":botPartCH, "buildOrientationNormal": botBuildOrientationNormal, "connectionArea": connectionAreaSum}) 
                     
                    # Removing unnecessary objects                
                    bpy.data.objects.remove(maxP["part"])     
                    bpy.data.objects.remove(maxP["CH"])              
                    convexParts.remove(maxP) 
                    manifoldCount =0
                    maxP=convexParts[0]  
                    for p in convexParts:
                        if max(p["part"].dimensions) > max(maxP["part"].dimensions):
                            maxP = p
                    planeLocation = maxP["part"].location   
                    planeNormal = plane_normal(maxP["part"])
                else:  
                    bpy.data.objects.remove(topPart) 
                    bpy.data.objects.remove(botPart) 
                    bpy.data.objects.remove(topPartCH) 
                    bpy.data.objects.remove(botPartCH) 
                    # Randomly moving plane location
                    rdLimit = maxP["part"].dimensions/10  
                    planeLocation = maxP["part"].location + Vector((rd.uniform(-rdLimit[0], rdLimit[0]), rd.uniform(-rdLimit[1], rdLimit[1]), rd.uniform(-rdLimit[2], rdLimit[2])))
                    manifoldCount += 1
                     
                bmTop.free()   
                bmBot.free()   
                 
                largeParts=False
                for p in convexParts:          
                    if max(p["part"].dimensions) > min(MAX_PART_SIZE):
                        largeParts=True
            decomposedParts.extend(convexParts)
        else:
            decomposedParts.append({"CH":i["CH"], "part":i["part"], "buildOrientationNormal": i["buildOrientationNormal"], "connectionArea": i["connectionArea"]})  
     
    # Checking a manifold geometry 
    print("===================================")
    for d in decomposedParts:
        manifoldAll = True
        bm = bmesh.new()
        bm.from_mesh(d["part"].data)   
        bm.edges.ensure_lookup_table()
        for e in bm.edges:
            if e.is_manifold==False:
                manifoldAll=False
        print(d["part"].name, " / ",manifoldAll, " / ",d["part"].dimensions)
        bm.free() 
            
    return decomposedParts, manifoldAll

# Post-processing of V-HACD: intersection of inputObj and its convex hull #######################################    
def intersect_inputObj_and_CH(input_obj, convexHulls):
    decomposedParts=[]
    manifoldAll=True
    for CH in convexHulls:     
        # Copy the original model for a boolean operation        
        bpy.ops.object.select_all(action='DESELECT')  
        input_obj.select_set(True)
        bpy.ops.object.duplicate()    
        copy_obj = bpy.context.selected_objects[0]    
        # Intersect operation between the original model and a convex hull
        b = copy_obj.modifiers.new(name='booly', type='BOOLEAN')
        b.object = CH
        b.operation = 'INTERSECT'
        bpy.ops.object.modifier_apply({"object": copy_obj},apply_as='DATA',modifier=b.name)     
          
        # Checking a manifold geometry
        bm = bmesh.new()
        bm.from_mesh(copy_obj.data)   
        bm.edges.ensure_lookup_table()
        for e in bm.edges:
            if e.is_manifold==False:
                manifoldAll=False
        bm.free() 
        if manifoldAll: 
            # Merging the faces of a convex hull with the same normal
            bmCH = bmesh.new()
            bmCH.from_mesh(CH.data) 
            bmCH.faces.ensure_lookup_table() 
            faces = bmCH.faces 
            faceMerge=[{"normal":faces[0].normal, "faceArea":faces[0].calc_area(), "index":[faces[0].index]}]
            for i in range(1, len(faces)):
                ck = False
                index = 0
                for j in range(0, len(faceMerge)):
                    if faceMerge[j]["normal"] == faces[i].normal:
                        ck = True
                        index = j
                if ck:
                    faceMerge[index]["faceArea"] += faces[i].calc_area()
                    faceMerge[index]["index"].append(faces[i].index)
                else:
                    faceMerge.append({"normal":faces[i].normal, "faceArea":faces[i].calc_area(), "index":[faces[i].index]})                 
            # Find the base plane (the face with the largest area)           
            largestArea=0
            for i in faceMerge: 
                if i["faceArea"] > largestArea:
                    largestArea=i["faceArea"] 
                    basePlaneNormal = i["normal"]   
                    largestFaceIndex=i["index"][0]
            mat = CH.matrix_world 
            tmp = mat @ -basePlaneNormal
            buildOrientationNormal = tmp.copy() 
            # Calculating the area of conection interface of a convex hull  
            connectionArea = 0 
            for f in bmCH.faces:  
                p = (f.calc_center_median()-copy_obj.location) * 0.999
                hit_data = copy_obj.ray_cast(p, f.normal)  
                if hit_data[0]: 
                    hitFaceIndex = hit_data[3]
                    bmPn = bmesh.new()
                    bmPn.from_mesh(copy_obj.data) 
                    bmPn.faces.ensure_lookup_table() 
                    connectionArea += bmPn.faces[hitFaceIndex].calc_area()  
                    bmPn.free()   
                      
            # Triangulating the decomposed parts
            bm = bmesh.new()
            bm.from_mesh(copy_obj.data) 
            bmesh.ops.triangulate(bm, faces=bm.faces[:], quad_method='BEAUTY', ngon_method='BEAUTY')
            bm.to_mesh(copy_obj.data)
            bm.free()
            bmCH.free()
            # Changing object color   
            bpy.context.view_layer.objects.active = copy_obj   
            bpy.context.object.color = (rd.random(), rd.random(), rd.random(), 1) 
            bpy.ops.object.origin_set( type = 'ORIGIN_GEOMETRY' )
            decomposedParts.append({"CH":CH, "part":copy_obj, "buildOrientationNormal": buildOrientationNormal, "connectionArea": connectionArea})
        else: 
            buildOrientationNormal = CH.matrix_world 
            decomposedParts.append({"CH":CH, "part":copy_obj, "buildOrientationNormal": buildOrientationNormal, "connectionArea": 0})
         
    return decomposedParts, manifoldAll


# Select an object   
input_obj = bpy.context.selected_objects[0]    
 
# Run concave decomposition based on V-HACD     
bpy.ops.object.vhacd() 

# Conserve the original shape of decomposed parts
convexHulls = bpy.context.selected_objects
decomposedParts, manifoldAll = intersect_inputObj_and_CH(input_obj, convexHulls)

# Run convex decomposition
if manifoldAll:
    decomposedParts, manifoldAll = PD_for_convex_features(decomposedParts)
else:
    print("A non-manifold geometry is detected") 

if manifoldAll:    
    # Obtain the evaluation indicators    
    P1 = roughness_of_part_surface(decomposedParts)
    P2 = intensity_of_critical_facets(decomposedParts)
    P3 = length_of_sharp_edges(decomposedParts)
    P4 = area_of_facets_causing_critical_gaps(decomposedParts)
    P5 = volume_of_concave_spaces(decomposedParts)
    f = feasibility_of_decomposed_parts(decomposedParts)
    if f == 0:
        P6 = "Feasible"
    else:
        P6 = "Unfeasible"        
    P7 = area_of_connection_interface(decomposedParts)
    P8 = len(decomposedParts)
     
    # Remove convex hulls
    for i in decomposedParts: 
        bpy.ops.object.select_all(action='DESELECT')  
        i["CH"].select_set(True)
        bpy.ops.object.delete()
        
    # Display the evaluation indicators    
    print("\n\n==========================================")
    print("Roughness: "+str(P1))
    print("Overhang: "+str(P2))
    print("Sharpness: "+str(P3)) 
    print("Gap: "+str(P4))
    print("Concavity: "+str(P5)) 
    print("Feasibility: "+str(P6))
    print("Interface: "+str(P7))
    print("Quantity: "+str(P8))
    print("==========================================")  
else:
    print("A non-manifold geometry is detected")

print("Computation time: ", timeit.default_timer() - tic)
