################################################################################################################### 
# Creator: Ha-An  
# Release: TBD
# Environment: Blender 2.82; Windows 10 Enterprise; 
# Requirement (add-on script): V-HACD (https://github.com/kmammou/v-hacd)
# Description: This code is a python script for Blender. This script represents the evaluation indicators of part decomposition for Additive Manufacturing. 
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
THRESHOLD_SHARP_EDGE = 170    # degrees
THRESHOLD_OVERHANG = 60             # degrees 
THRESHOLD_GAP_DISTANCE = 2          # meters : 1
MAX_PART_SIZE = [200, 200, 300] # meters 
THRESHOLD_MIN_SIZE = 1      # meters 
# Computation time
tic = timeit.default_timer() 


'''
def is_watertight(object: bpy.types.Object, check_self_intersection=True) -> bool:
    """
    Checks whether the given object is watertight or not
    :param object: Object the inspect
    :return: True if watertight, False otherwise
    """
    old_active_object = bpy.context.view_layer.objects.active
    old_mode = old_active_object.mode
    bpy.context.view_layer.objects.active = object
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_non_manifold(extend=False)
    bm = bmesh.from_edit_mesh(object.data
    is_watertight = True
    for v in bm.verts:
        if v.select:
            is_watertight = False
            break
    # Only check for self intersection if manifold
    if is_watertight and check_self_intersection:
        bvh_tree = mathutils.bvhtree.BVHTree.FromBMesh(bm, epsilon=0.000001)
        intersections = bvh_tree.overlap(bvh_tree)
        if intersections:
            is_watertight = False
    bpy.context.view_layer.objects.active = old_active_object
    bpy.ops.object.mode_set(mode=old_mode)
    return is_watertight
'''
# Demir (2018) properties ############################################################################################   
def concavity(obj, CH): 
    bpy.ops.object.select_all(action='DESELECT')   
    CH.select_set(True) 
    bpy.ops.object.origin_set( type = 'ORIGIN_GEOMETRY' )
    concavityMeasurement=0
    # Creating mesh data for the input object (obj) 
    bmObj = bmesh.new()
    bmObj.from_mesh(obj.data)
    bmObj.verts.ensure_lookup_table() 
    # Creating mesh data for a convex hull (CH)
    bmCH = bmesh.new()
    bmCH.from_mesh(CH.data)
    bmCH.faces.ensure_lookup_table() 
    # Switching matrix
    mat_global_obj = obj.matrix_world # Local (obj) to global
    mat_local_obj = obj.matrix_world.inverted() # Global to local (obj)
    mat_global_CH = CH.matrix_world   # Local (CH) to global
    mat_local_CH = CH.matrix_world.inverted()   # Global to local (CH)  
#    VV=[Vector((0,0,0)), 0, 0]
#    print("obj: ", obj.name)
    for v in bmObj.verts: # for all vertex of the original model
        global_obj_v = mat_global_obj @ v.co # a vertex of an object (in global space)
        local_CH_v = mat_local_CH @ global_obj_v # a vertex of an object (in local space for the convex hull)     
        for f in bmCH.faces:   
            hitCH_data = CH.ray_cast(local_CH_v, f.normal) # hitData structure: (result, location, normal, face index)  
            if hitCH_data[0]: # Did the ray hit the face of CH?  
                if hitCH_data[3] == f.index and hitCH_data[1] != local_CH_v:   
                    distance_obj_to_CH = distance.euclidean(hitCH_data[1] , local_CH_v) # Project the length of delta along normal 
                    
                    global_hitCH_point = mat_global_CH @ (hitCH_data[1] + f.normal*0.1 )
                    global_hitCH_normal = mat_global_CH @ hitCH_data[2]   
                    local_obj_v = mat_local_obj @ global_hitCH_point # Converting global space into object space (Object) 
                    local_obj_normal = mat_local_obj @ global_hitCH_normal # Converting global space into object space (Object)  
                     
                    hitObj_data = obj.ray_cast(local_obj_v, -local_obj_normal) # hitData structure: (result, location, normal, face index)          
                    if hitObj_data[0] and hitObj_data[1] != local_obj_v:  
                        distance_CH_to_obj = distance.euclidean(hitObj_data[1], local_obj_v) # Project the length of delta along normal.    
                        
                        if abs(distance_obj_to_CH - distance_CH_to_obj) < 1 :  
                            if distance_CH_to_obj>concavityMeasurement: 
                                concavityMeasurement=distance_obj_to_CH 
#                                largestFace = f.copy()   
#                                VV[0] = Vector((0,0,0))
#                                VV[1] = hitCH_data[1].copy()
#                                VV[2] = local_CH_v.copy()
    '''  
    # Leave non-allowable faces as objects for validation
    mesh = bpy.data.meshes.new("f")  # add the new mesh
    obj_f = bpy.data.objects.new(mesh.name, mesh)
    bpy.context.collection.objects.link(obj_f)
    bpy.context.view_layer.objects.active = obj_f
    v = [] 
    f = [[]]
    mat_global_obj = CH.matrix_world
    for i in range(0, len(largestFace.verts)):  
        v.append(mat_global_obj @ largestFace.verts[i].co)
        f[0].append(i)
    edges = []
    mesh.from_pydata(v, [], f) 
    # Leave non-allowable faces as objects for validation
    mesh = bpy.data.meshes.new("d")  # add the new mesh
    obj_f = bpy.data.objects.new(mesh.name, mesh)
    bpy.context.collection.objects.link(obj_f)
    bpy.context.view_layer.objects.active = obj_f
    v = [] 
    f = [[]]
    mat_global_obj = CH.matrix_world
    for i in range(0, len(VV)):  
        v.append(mat_global_obj @ VV[i])
        f[0].append(i)
    edges = []
    mesh.from_pydata(v, [], f) 
    ''' 
#    print("concavityMeasurement ", concavityMeasurement)
    return concavityMeasurement

def surface_angle(obj, buildOrientationNormal):
    surfaceAngleSum =0 
    bmObj = bmesh.new()
    bmObj.from_mesh(obj.data)
    bmObj.faces.ensure_lookup_table() 
    # Calculate surface angles
    for f in bmObj.faces:   
        if f.normal != Vector((0.0, 0.0, 0.0)):
            a = degrees(buildOrientationNormal.angle(f.normal)) 
            if a<=THRESHOLD_SURFACE_ROUGHNESS:
                surfaceAngleSum += a  
                '''
                # Leave non-allowable faces as objects for validation
                mesh = bpy.data.meshes.new("m")  # add the new mesh
                obj_f = bpy.data.objects.new(mesh.name, mesh)
                bpy.context.collection.objects.link(obj_f)
                bpy.context.view_layer.objects.active = obj_f
                v = [] 
                faces = [[]]
                mat_global_obj = obj.matrix_world
                for i in range(0, len(f.verts)):  
                    v.append(mat_global_obj @ f.verts[i].co)
                    faces[0].append(i)
                edges = []
                mesh.from_pydata(v, [], faces) 
                '''
    return surfaceAngleSum
    
# The first property for concavity
def concavity_property(decomposedParts): 
    concavityList=[]    
    for d in decomposedParts:        
        concavityList.append(concavity(d["part"], d["CH"]))    
    return mean(concavityList)
# The second property for surface angles
def surface_angle_property(decomposedParts):    
    surfaceAngle=[]    
    for d in decomposedParts: 
        surfaceAngle.append(surface_angle(d["part"], d["buildOrientationNormal"])) 
    return mean(surfaceAngle)
# Third property for the balance between size and number of components
def size_balance(decomposedParts):
    sizeBalance=[]
    list_dimensionX=[]
    list_dimensionY=[]
    list_dimensionZ=[]
    for d in decomposedParts: 
#        print(d["part"].name, ": ", d["part"].dimensions) 
        list_dimensionX.append(d["part"].dimensions.x)
        list_dimensionY.append(d["part"].dimensions.y)
        list_dimensionZ.append(d["part"].dimensions.z)
    meanDimension=[mean(list_dimensionX), mean(list_dimensionY), mean(list_dimensionZ)]
    for d in decomposedParts:  
        deviationX = d["part"].dimensions.x - meanDimension[0]
        deviationY = d["part"].dimensions.y - meanDimension[1]
        deviationZ = d["part"].dimensions.z - meanDimension[2]
        p1 = (0, 0, 0)
        p2 = (deviationX, deviationY, deviationZ) 
        sizeBalance.append(distance.euclidean(p1, p2)) # Calculate euclidean distance         
    return mean(sizeBalance)
# The fourth property for shape deviation (based-based)
def shape_deviation(decomposedParts):
    shapeVolumeList=[]
    for d in decomposedParts:
        part = d["part"]
        CH = d["CH"]    
        bmPart = bmesh.new()
        bmPart.from_mesh(part.data)   
        volumePart = bmPart.calc_volume(signed=True)  
        bmCH = bmesh.new()
        bmCH.from_mesh(CH.data)   
        volumeCH = bmCH.calc_volume(signed=True)   
#        print("part:", d["part"].name, ", Vol: ",max((volumeCH - volumePart), 0))
        shapeVolumeList.append(max((volumeCH - volumePart), 0))  
    return mean(shapeVolumeList)  
# ROUHGNESS #####################################################################################################
def roughness_for_each_part(part, CH, buildOrientationNormal):     
    # Calculate surface roughness
    bmPart = bmesh.new()
    bmPart.from_mesh(part.data) 
    surfaceRoughness =0
    bmPart.faces.ensure_lookup_table() 
    for f in bmPart.faces:  
#        mat = part.matrix_world 
#        global_f_normal = mat @ f.normal
        global_f_normal = f.normal
        if global_f_normal != Vector((0.0, 0.0, 0.0)):
            theta = buildOrientationNormal.angle(global_f_normal)   
#            print("f.index: ",f.index) 
#            print("theta: ",degrees(theta))    
            if 0<theta<radians(180):   
    #            print("surfaceAngle: ", degrees(theta))
    #            print("surfaceArea: ", f.calc_area())
                SR = abs(cos(theta))
            else: 
                SR = 0
    #        print("surfaceAngle: ", degrees(theta))
    #        print("SR: ",SR) 
    #        print("Area: ", f.calc_area())
            surfaceRoughness += SR * f.calc_area()    
    bmPart.free()  
#    print("surfaceRoughness: ", surfaceRoughness)
    return surfaceRoughness  
def roughness_of_part_surface(decomposedParts): 
    sum=0    
    for i in decomposedParts:    
#        print("Part: ", i["part"].name, "; CH: ", i["CH"].name)    
        sum += roughness_for_each_part(i["part"], i["CH"], i["buildOrientationNormal"])
    return sum 

# INTENSITY   ##########################################################################################
def intensity_of_critical_facets_for_each_part(part, CH, buildOrientationNormal):  
    # Calculate overhang area
    bmPart = bmesh.new()
    bmPart.from_mesh(part.data) 
    unsupportedSurfaceIntensity =0 
    bmPart.faces.ensure_lookup_table()  
#    print("buildOrientationNormal: ", buildOrientationNormal)
    for f in bmPart.faces:  
#        mat = part.matrix_world 
#        global_f_normal = mat @ f.normal
        global_f_normal = f.normal
        if global_f_normal != Vector((0.0, 0.0, 0.0)):
            theta = buildOrientationNormal.angle(global_f_normal)   
            if radians(90) < theta:
                if radians(180 - THRESHOLD_OVERHANG) < theta and theta < radians(180-0.05):  ## 0.5 is an adjustment value for implementation
#                    print("surfaceAngle: ", degrees(theta))
        #            print("surfaceArea: ", f.calc_area())
                    OH = abs(cos(theta))   
                    '''
                    # Leave non-allowable faces as objects for validation
                    mesh = bpy.data.meshes.new("m")  # add the new mesh
                    obj_f = bpy.data.objects.new(mesh.name, mesh)
                    bpy.context.collection.objects.link(obj_f)
                    bpy.context.view_layer.objects.active = obj_f
                    verts = [] 
                    faces = [[]]  
                    mat = part.matrix_world   
                    for i in range(0, len(f.verts)): 
                        verts.append(mat @ f.verts[i].co)
                        faces[0].append(i)
                    edges = []
                    mesh.from_pydata(verts, [], faces)
                    '''
                else: 
                    OH = 0
        #        print("surfaceAngle: ", degrees(theta))
        #        print("SR: ",SR) 
        #        print("Area: ", f.calc_area())
                unsupportedSurfaceIntensity += OH * f.calc_area()     
    bmPart.free() 
#    print("unsupportedSurfaceIntensity: ", unsupportedSurfaceIntensity)
    return unsupportedSurfaceIntensity
# Sum of the area of causing overhangs for all convex hulls
def intensity_of_critical_facets(decomposedParts): 
    sum=0    
    for i in decomposedParts:    
#        print("Part: ", i["part"].name, "; CH: ", i["CH"].name)   
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
#            print("theta: ", degrees(f1.normal.angle(f2.normal)))
            if theta > THRESHOLD_SHARP_EDGE:
                sharpEdgeLength += e.calc_length()  
                '''
                # Leave non-allowable faces as objects for validation
                mesh = bpy.data.meshes.new("m")  # add the new mesh
                obj_f = bpy.data.objects.new(mesh.name, mesh)
                bpy.context.collection.objects.link(obj_f)
                bpy.context.view_layer.objects.active = obj_f
                verts = [] 
                faces = [[]]  
                mat = part.matrix_world   
                for i in range(0, len(f1.verts)): 
                    verts.append(mat @ f1.verts[i].co)
                    faces[0].append(i)
                edges = []
                mesh.from_pydata(verts, [], faces)
                # Leave non-allowable faces as objects for validation
                mesh = bpy.data.meshes.new("m")  # add the new mesh
                obj_f = bpy.data.objects.new(mesh.name, mesh)
                bpy.context.collection.objects.link(obj_f)
                bpy.context.view_layer.objects.active = obj_f
                verts = [] 
                faces = [[]]  
                mat = part.matrix_world   
                for i in range(0, len(f2.verts)): 
                    verts.append(mat @ f2.verts[i].co)
                    faces[0].append(i)
                edges = []
                mesh.from_pydata(verts, [], faces)
                '''
                
    bmPart.free()  
#    print("surfaceRoughness: ", surfaceRoughness)
    return sharpEdgeLength  
def length_of_sharp_edges(decomposedParts): 
    sum=0    
    for i in decomposedParts:    
#        print("Part: ", i["part"].name, "; CH: ", i["CH"].name)    
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
#                print(f.calc_area())
                nonAllowableGapArea +=f.calc_area()    
                '''
                # Leave non-allowable faces as objects for validation
                mesh = bpy.data.meshes.new("m")  # add the new mesh
                obj_f = bpy.data.objects.new(mesh.name, mesh)
                bpy.context.collection.objects.link(obj_f)
                bpy.context.view_layer.objects.active = obj_f
                v = [] 
                faces = [[]]
                mat_global_obj = part.matrix_world
                for i in range(0, len(f.verts)):  
                    v.append(mat_global_obj @ f.verts[i].co)
                    faces[0].append(i)
                edges = []
                mesh.from_pydata(v, [], faces)  
                '''
                
    '''
    for i in range(0, len(bmObj.faces)):  
        # distance check by using RayCast
        f = bmObj.faces[i]   
        p = f.calc_center_median()+f.normal*0.001 
        hit_data = part.ray_cast(p, f.normal)     
        gapDistance=0
        if hit_data[0]:    
            gapDistance = distance.euclidean(f.calc_center_median(), hit_data[1]) 
        if gapDistance>0:
            if gapDistance < THRESHOLD_GAP_DISTANCE:
                bpy.ops.object.select_all(action='DESELECT') # Deselect all objects
#                print(f.calc_area())
                nonAllowableGapArea +=f.calc_area()        
    '''    
    bmObj.free()                
#    print("nonAllowableGapArea: ", nonAllowableGapArea)    
    return nonAllowableGapArea 
def area_of_facets_causing_critical_gaps(decomposedParts): 
    sum=0    
    for i in decomposedParts:    
#        print("Part: ", i["part"].name, "; CH: ", i["CH"].name)   
        a = area_of_facets_causing_critical_gaps_for_each_part(i["part"])
        sum += a
    return sum 

'''
# Calculate the curvature of the critical surface ############################################## 
def surface_curvature(part): 
    bm = bmesh.new()
    bm.from_mesh(part.data) 
    bm.verts.ensure_lookup_table()   
    bm.faces.ensure_lookup_table()       
    criticalCurvature=0 
    for v in bm.verts: 
        print(v)
        # Calculating AV 
        faces=v.link_faces 
        sum=0
        for f in faces: 
            sum += f.calc_area() 
        AV = sum/3 
        # Calculating the sum of theta 
        angle=0 #the sum of theta 
        for f in faces:
#            print(f) 
            sv=[]   
            for i in f.verts: 
                if i.index != v.index:
        #            print(i)
                    sv.append(i)
            v1 = sv[0].co - v.co
            v2 = sv[1].co - v.co
            v1mag = math.sqrt(v1.x * v1.x + v1.y * v1.y + v1.z * v1.z)
            v1norm = [v1.x/v1mag , v1.y/v1mag , v1.z/v1mag]
            v2mag = math.sqrt(v2.x * v2.x + v2.y * v2.y + v2.z * v2.z)
            v2norm = [v2.x/v2mag , v2.y/v2mag , v2.z/v2mag]
            res = v1norm[0] * v2norm[0] + v1norm[1] * v2norm[1] + v1norm[2] * v2norm[2]
            theta = (math.acos(res)*180/math.pi)     
#            print(theta)
            angle += radians(theta)
        # Calculating the Gaussian curvature 
        #print(degrees(angle)) 
        #GC = (2*radians(180)-angle) / AV
        GC = angle / (2*radians(180))
        print(round(GC, 2))
        if GC < THRESHOLD_CURVATURE:
            GC = 0  
        criticalCurvature += GC 
    bm.free()   
    return criticalCurvature
# Sum of surface curvature
def sum_surface_curvature(decomposedParts): 
    sum=0    
    for i in decomposedParts:        
#        print("Part: ", i["part"].name, "; CH: ", i["CH"].name)
        v=surface_curvature(i["part"])
#        print(v)
        sum += v
    return sum
''' 
# CONCAVITY #########################################################################################################
def volume_of_concave_spaces_for_each_part(part, CH): 
    bmPart = bmesh.new()
    bmPart.from_mesh(part.data)   
    volumePart = bmPart.calc_volume(signed=True)  
    bmCH = bmesh.new()
    bmCH.from_mesh(CH.data)   
    volumeCH = bmCH.calc_volume(signed=True)   
#    print("volumeCH: ", volumeCH)
#    print("volumePart: ", volumePart)
    return max((volumeCH - volumePart), 0) 
# Sum of concave space volume for all convex hulls
def volume_of_concave_spaces(decomposedParts): 
    sum=0    
    for i in decomposedParts:        
#        print("Part: ", i["part"].name, "; CH: ", i["CH"].name)
        v=volume_of_concave_spaces_for_each_part(i["part"], i["CH"])
#        print(v)
        sum += v
    return sum 
 
# INTERFACE #################################################################################################
'''
def contact_area(part, CH):  
    mat_global_obj = part.matrix_world  
    # Creating mesh data for the input part
    bmCH = bmesh.new()
    bmCH.from_mesh(CH.data)
    bmCH.faces.ensure_lookup_table()    
    cuttingAreaSum = 0
    for f in bmCH.faces:  
        p = f.calc_center_median()-f.normal*0.001 
        hit_data = part.ray_cast(p, f.normal)  
        if hit_data[0]: 
            # Finding contact faces
            mesh = bpy.data.meshes.new("m")  # add the new mesh
            obj_f = bpy.data.objects.new(mesh.name, mesh)
            bpy.context.collection.objects.link(obj_f)
            bpy.context.view_layer.objects.active = obj_f 
            verts = [] 
            faces =[[]]  
            for i in range(0, len(f.verts)): 
                verts.append(mat_global_obj @ f.verts[i].co)  
                faces[0].append(i)
            edges = [] 
            mesh.from_pydata(verts, [], faces)  
            # Projection by using Shrinkwrap 
            bpy.ops.object.select_all(action='DESELECT')  
            part.select_set(True)
            bpy.ops.object.duplicate()    
            copy_obj = bpy.context.selected_objects[0]   
            bpy.ops.object.select_all(action='DESELECT')   
            copy_obj.select_set(True)  
            bpy.context.view_layer.objects.active = copy_obj
            bpy.ops.object.modifier_add(type='SHRINKWRAP')
            bpy.context.object.modifiers["Shrinkwrap"].target = obj_f
            bpy.context.object.modifiers["Shrinkwrap"].wrap_method = 'TARGET_PROJECT'
            bpy.ops.object.modifier_apply(apply_as='DATA', modifier="Shrinkwrap")
                 
            bmA = bmesh.new()
            bmA.from_mesh(copy_obj.data) 
            bmA.faces.ensure_lookup_table()    
            cuttingAreaSum += sum([f.calc_area() for f in bmA.faces])
            
            bpy.data.objects.remove(obj_f)  
            bpy.data.objects.remove(copy_obj)  
            bmA.free()
    bmCH.free()          
    return cuttingAreaSum 
''' 
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
        #global_dim = i["part"].matrix_world @ i["part"].dimensions 
        maxPartDim = max(i["part"].dimensions)
#        maxPartDim = math.sqrt(i["part"].dimensions.x**2+i["part"].dimensions.y**2+i["part"].dimensions.z**2)
        minWorkspaceDim = min(MAX_PART_SIZE)  
#        print("i: ", i["part"])   
#        print("maxPartDim: ", maxPartDim)      
#        print("minWorkspaceDim: ", minWorkspaceDim)      
        if maxPartDim < minWorkspaceDim:
            M = 0
        else:
            M = inf 
            print("i: ", i["part"])   
            print("maxPartDim: ", maxPartDim)      
            print("minWorkspaceDim: ", minWorkspaceDim)    
            break
        '''
        minPartDim = min(i["part"].dimensions)
        minBuildableDim = THRESHOLD_MIN_SIZE 
        if minPartDim > minBuildableDim:
            M = 0
        else:
            M = inf 
            print("i: ", i["part"])   
            print("minPartDim: ", minPartDim)      
            print("minBuildableDim: ", minBuildableDim)    
            break
        '''
    return M 
'''
# Calculating size deviation among decomposed parts ############################################################
def size_deviation_among_parts_3(decomposedParts): 
    listDimensionX=[]
    listDimensionY=[]
    listDimensionZ=[]
    for i in decomposedParts:     
        global_dim = i["part"].matrix_world @ i["part"].dimensions
        listDimensionX.append(global_dim.x)
        listDimensionY.append(global_dim.y)
        listDimensionZ.append(global_dim.z)
    meanDimension=[mean(listDimensionX), mean(listDimensionY), mean(listDimensionZ)]
    sizeDeviation=[]
    for i in decomposedParts:    
        global_dim = i["part"].matrix_world @ i["part"].dimensions
        deviationX = global_dim.x - meanDimension[0]
        deviationY = global_dim.y - meanDimension[1]
        deviationZ = global_dim.z - meanDimension[2]
        p1 = (0, 0, 0)
        p2 = (deviationX, deviationY, deviationZ) 
#        print("distance.euclidean(p1, p2): ", distance.euclidean(p1, p2))
        sizeDeviation.append(distance.euclidean(p1, p2)) # Calculate euclidean distance  
#    print("mean(sizeDeviation): ", mean(sizeDeviation))    
    return mean(sizeDeviation)
''' 
# Part decomposition for convex featuresn ######################################################################################
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
#    print("buildOrientationNormal ", buildOrientationNormal)
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
                
#                for c in convexParts:
#                    print("convexParts 1", c["part"].name)
#                    print("dimensions ", c["part"].dimensions) 
                print("maxP ", maxP["part"].name) 
                print("dimensions ", maxP["part"].dimensions)  
   
####################bpy.context.active_object.dimensions = maxP["part"].dimensions  ##########VALIDATION#########
                       
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
#                    print("convexParts 2", p["part"].name) 
#                    print("dimensions ", p["part"].dimensions)         
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
            ''' 
            # Leave non-allowable faces as objects for validation
            mesh = bpy.data.meshes.new("base")  # add the new mesh
            obj_f = bpy.data.objects.new(mesh.name, mesh)
            bpy.context.collection.objects.link(obj_f)
            bpy.context.view_layer.objects.active = obj_f
            v = [] 
            f = [[]]
            mat_global_obj = CH.matrix_world
            for i in range(0, len(faces[largestFaceIndex].verts)):  
                v.append(mat_global_obj @ faces[largestFaceIndex].verts[i].co)
                f[0].append(i)
            edges = []
            mesh.from_pydata(v, [], f) 
            '''
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
 
# Run V-HACD   0.05106692 0.0253298  0.         0.       
bpy.ops.object.vhacd()
#bpy.ops.object.vhacd(gamma=1.0)     
#bpy.ops.object.vhacd(resolution=50000) 
#bpy.ops.object.vhacd( alpha=0.05106692, beta=0.0253298, gamma=0, concavity= 0)
#bpy.ops.object.vhacd(resolution=55.56608552*10000, minVolumePerCH=21.2288925*0.0001, maxNumVerticesPerCH=1.49137808*4) 
# Run post-processing 
convexHulls = bpy.context.selected_objects
decomposedParts, manifoldAll = intersect_inputObj_and_CH(input_obj, convexHulls)

'''
# Run convex decomposition
if manifoldAll:
    decomposedParts, manifoldAll = PD_for_convex_features(decomposedParts)
else:
    print("A non-manifold geometry is detected")
'''
'''
# Export the decomposed models 
path='C://Users/jspoh/Google Drive/03.Research/8.Journal/2020AM (PD)/Data/CaseStudy/GoldenGateBridge/Demir(2018)/'
for i in decomposedParts:
    bpy.ops.object.select_all(action='DESELECT') 
    i["part"].select_set(True)
    bpy.ops.export_mesh.stl(filepath=path+i["part"].name+".stl", use_selection=True)   
    bpy.ops.object.select_all(action='DESELECT') 
    i["CH"].select_set(True)
    bpy.ops.export_mesh.stl(filepath=path+i["CH"].name+".stl", use_selection=True)   
'''
'''
input_obj = bpy.context.scene.objects["SB low"] 
convexHulls = [bpy.context.scene.objects["SB low_hull_1"]]
decomposedParts=intersect_inputObj_and_CH(input_obj, convexHulls)
#decomposedParts.append({"part":bpy.context.scene.objects["Icosahedron"], "CH":bpy.context.scene.objects["Icosahedron.001"]})
'''

# Changing object color 
for i in decomposedParts:
    bpy.ops.object.select_all(action='DESELECT')  
    i["CH"].select_set(True)
    bpy.context.view_layer.objects.active = i["CH"]  
    bpy.context.object.display_type = 'TEXTURED'
    bpy.context.object.color = (rd.random(), rd.random(), rd.random(), 1) 


'''
if manifoldAll: 
    P1 = roughness_of_part_surface(decomposedParts)
    P2 = intensity_of_critical_facets(decomposedParts)
    P3 = length_of_sharp_edges(decomposedParts)
    P4 = area_of_facets_causing_critical_gaps(decomposedParts)
    P5 = volume_of_concave_spaces(decomposedParts)
    P6 = area_of_connection_interface(decomposedParts)
    P7 = feasibility_of_decomposed_parts(decomposedParts)
    P8 = len(decomposedParts)
    print("==========================================")
    print("roughness_of_part_surface: "+str(P1))
    print("intensity_of_critical_facets: "+str(P2))
    print("length_of_sharp_edges: "+str(P3)) 
    print("area_of_facets_causing_critical_gaps: "+str(P4))
    print("volume_of_concave_spaces: "+str(P5)) 
    print("area_of_connection_interface: "+str(P6))
    print("feasibility_of_decomposed_parts: "+str(P7))
    print("number_of_decomposed_parts: "+str(P8))
    print("==========================================") 
    print("==========================================")
    print(P1)
    print(P2)
    print(P3)
    print(P4)
    print(P5)
    print(P6)
    print(P7)
    print(P8)
    print("==========================================") 
    
    
    # The previous work of Demir (2018) ###################################### 
    print("==========================================") 
    print("concavity_property:  ", concavity_property(decomposedParts))
    print("surface_angle_property:  ", surface_angle_property(decomposedParts))
    print("size_balance:  ", size_balance(decomposedParts))
    print("shape_deviation:  ", shape_deviation(decomposedParts) ) 
    print("==========================================")
    
else:
    print("A non-manifold geometry is detected")
'''
'''
# File validation
input_obj = bpy.context.selected_objects[0]  
val = input_obj.data.validate() # Validate geometry: return True when the mesh has had invalid geometry corrected/removed
wat = is_watertight(input_obj) # Validate watertight 
print("VAL-", val) 
print("WAT-", wat)
'''

print("Computation time: ", timeit.default_timer() - tic)
