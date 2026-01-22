"""
McStas Instrument Visualizer for Jupyter Notebooks using pythreejs
"""

import re
import numpy as np
import pythreejs as p3
import ipywidgets as ipw
import math


SKIP_TYPES = [
    "arm",
    "union_init",
    "union_master",
    "union_stop",
    "union_make_material",
    "incoherent_process",
    "powder_process",
]


class McStasComponent:
    """Represents a single McStas component with position and rotation"""

    def __init__(
        self,
        name,
        component_type,
        at_pos,
        relative_to="ABSOLUTE",
        rotation=None,
        params=None,
    ):
        self.name = name
        self.component_type = component_type
        self.at_pos = at_pos  # [x, y, z]
        self.relative_to = relative_to
        self.rotation = (
            rotation if rotation is not None else [0, 0, 0]
        )  # [rx, ry, rz] in degrees
        self.params = params if params is not None else {}
        self.absolute_pos = None
        self.otation = None

    def __repr__(self):
        return f"McStasComponent({self.name}, {self.component_type}, pos={self.at_pos}, rel={self.relative_to})"


class McStasInstrumentParser:
    """Parses McStas instrument files"""

    def __init__(self, instrument_file):
        self.instrument_file = instrument_file
        self.components = []
        self.instrument_name = ""
        self.parameters = {}

    def parse(self):
        """Parse the instrument file"""
        with open(self.instrument_file, "r") as f:
            content = f.read()

        # Extract instrument name
        instrument_match = re.search(r"DEFINE\s+INSTRUMENT\s+(\w+)", content)
        if instrument_match:
            self.instrument_name = instrument_match.group(1)

        # # Extract and parse instrument parameters with default values
        # params_match = re.search(
        #     r"DEFINE\s+INSTRUMENT\s+\w+\s*\((.*?)\)", content, re.DOTALL
        # )
        # if params_match:
        #     params_str = params_match.group(1)
        #     # Parse parameter = value pairs
        #     param_matches = re.finditer(r"double\s+(\w+)\s*=\s*([^,\)]+)", params_str)
        #     for match in param_matches:
        #         param_name = match.group(1).strip()
        #         param_value = match.group(2).strip()
        #         try:
        #             self.parameters[param_name] = float(param_value)
        #         except ValueError:
        #             self.parameters[param_name] = 0.0  # Default if can't parse

        # # Extract and parse instrument parameters with default values
        # params_match = re.search(
        #     r"DEFINE\s+INSTRUMENT\s+\w+\s*\((.*?)\)", content, re.DOTALL
        # )
        # if params_match:
        #     params_str = params_match.group(1)
        #     print(f"DEBUG PARAMS STRING:\n{params_str}\n{'=' * 80}")
        #     # Parse parameter = value pairs
        #     param_matches = re.finditer(r"double\s+(\w+)\s*=\s*([^,\)]+)", params_str)
        #     for match in param_matches:
        #         param_name = match.group(1).strip()
        #         param_value = match.group(2).strip()
        #         try:
        #             self.parameters[param_name] = float(param_value)
        #         except ValueError:
        #             self.parameters[param_name] = 0.0  # Default if can't parse
        #         print(
        #             f"DEBUG PARAM:  {param_name} = {param_value} -> {self.parameters[param_name]}"
        #         )

        # print(f"DEBUG:  Instrument parameters: {self.parameters}")

        # Extract and parse instrument parameters with default values
        params_match = re.search(
            r"DEFINE\s+INSTRUMENT\s+\w+\s*\((.*?)\)", content, re.DOTALL
        )
        if params_match:
            params_str = params_match.group(1)
            print(f"DEBUG PARAMS STRING:\n{params_str}\n{'=' * 80}")
            # Parse parameter = value pairs
            param_matches = re.finditer(r"double\s+(\w+)\s*=\s*([^,\)]+)", params_str)
            for match in param_matches:
                param_name = match.group(1).strip()
                param_value = match.group(2).strip()
                # Remove comments (// style)
                param_value = re.sub(r"//.*$", "", param_value).strip()
                try:
                    self.parameters[param_name] = float(param_value)
                except ValueError:
                    self.parameters[param_name] = 0.0  # Default if can't parse
                print(
                    f"DEBUG PARAM:  {param_name} = {param_value} -> {self.parameters[param_name]}"
                )

        print(f"DEBUG:  Instrument parameters: {self.parameters}")

        # Extract and parse DECLARE section variables
        declare_match = re.search(r"DECLARE\s*%\{(.*?)%\}", content, re.DOTALL)
        if declare_match:
            declare_str = declare_match.group(1)
            print(f"DEBUG DECLARE STRING:\n{declare_str[:500]}\n{'=' * 80}")
            # Parse double variable = value pairs
            declare_matches = re.finditer(r"double\s+(\w+)\s*=\s*([^;]+);", declare_str)
            for match in declare_matches:
                var_name = match.group(1).strip()
                var_value = match.group(2).strip()
                # Remove comments (// style)
                var_value = re.sub(r"//.*$", "", var_value).strip()
                try:
                    self.parameters[var_name] = float(var_value)
                    print(
                        f"DEBUG DECLARE VAR: {var_name} = {var_value} -> {self.parameters[var_name]}"
                    )
                except ValueError:
                    # Might be an expression or string, skip for now
                    print(
                        f"DEBUG DECLARE VAR: {var_name} = {var_value} (could not parse as float)"
                    )
                    pass

        print(f"DEBUG: All parameters after DECLARE:  {self.parameters}")

        # print(f"DEBUG: Instrument parameters: {self.parameters}")

        # Extract TRACE section
        trace_match = re.search(
            r"TRACE\s+(.*?)\s+(?:FINALLY|SAVE|END)\s",
            content,
            re.DOTALL | re.IGNORECASE,
        )
        if not trace_match:
            raise ValueError("Could not find TRACE section in instrument file")

        trace_section = trace_match.group(1)

        # print(f"DEBUG: Full TRACE section:\n{trace_section}\n{'=' * 80}")

        # Parse components - updated pattern to handle multi-line parameters and EXTEND blocks
        # Split by COMPONENT keyword first
        component_blocks = re.split(r"\bCOMPONENT\s+", trace_section)[
            1:
        ]  # Skip first empty split

        # print(f"DEBUG: Found {len(component_blocks)} component blocks after split")

        for i, block in enumerate(component_blocks):
            # print(f"\nDEBUG: Block {i}:\n{block[:200]}.. .\n{'-' * 40}")

            # Extract component name and type - match up to AT keyword
            header_match = re.match(
                r"(\w+)\s*=\s*(\w+)\s*\((.*?)\)\s*", block, re.DOTALL
            )

            # if not header_match:
            #     print(f"DEBUG: Block {i} - header_match FAILED")
            #     continue

            name = header_match.group(1)
            comp_type = header_match.group(2)
            params_str = header_match.group(3)

            # Find the AT clause (may come after WHEN or directly)
            at_match = re.search(
                r"AT\s*\((.*?)\)\s*(RELATIVE|ABSOLUTE)\s*(\w*)", block, re.DOTALL
            )

            # if not at_match:
            #     print(f"DEBUG: Block {i} ({name}) - at_match FAILED")
            #     continue

            at_str = at_match.group(1)
            rel_type = at_match.group(2)
            rel_to = (
                at_match.group(3).strip() if at_match.group(3).strip() else "ABSOLUTE"
            )

            # print(
            #     f"DEBUG: Found component:   {name}, type: {comp_type}, at: {at_str.strip()}, rel: {rel_to}"
            # )

            # # Parse position - evaluate expressions using instrument parameters
            # at_pos = []
            # for coord in at_str.split(","):
            #     coord = coord.strip()
            #     try:
            #         # Try direct conversion first
            #         at_pos.append(float(coord))
            #     except ValueError:
            #         # Try to evaluate as expression with instrument parameters
            #         try:
            #             # Create evaluation context with parameters and math functions
            #             import math

            #             eval_context = {
            #                 **self.parameters,
            #                 "PI": math.pi,
            #                 "sin": math.sin,
            #                 "cos": math.cos,
            #                 "tan": math.tan,
            #                 "atan": math.atan,
            #             }
            #             value = eval(coord, {"__builtins__": {}}, eval_context)
            #             at_pos.append(float(value))
            #         except:
            #             print(f"WARNING: Could not evaluate '{coord}', using 0.0")
            #             at_pos.append(0.0)

            eval_context = {
                **self.parameters,
                "PI": math.pi,
                "sin": math.sin,
                "cos": math.cos,
                "tan": math.tan,
                "atan": math.atan,
            }

            # Parse position - evaluate expressions using instrument parameters
            at_pos = []
            for coord in at_str.split(","):
                coord = coord.strip()
                try:
                    # Try direct conversion first
                    at_pos.append(float(coord))
                except ValueError:
                    # Try to evaluate as expression with instrument parameters
                    try:
                        # # Create evaluation context with parameters and math functions
                        # import math

                        # eval_context = {
                        #     **self.parameters,
                        #     "PI": math.pi,
                        #     "sin": math.sin,
                        #     "cos": math.cos,
                        #     "tan": math.tan,
                        #     "atan": math.atan,
                        # }
                        value = eval(coord, {"__builtins__": {}}, eval_context)
                        at_pos.append(float(value))
                        print(f"DEBUG EVAL: {name} - evaluated '{coord}' = {value}")
                    except Exception as e:
                        print(
                            f"WARNING: Could not evaluate '{coord}' for {name}: {e}, using 0.0"
                        )
                        at_pos.append(0.0)

            print(f"DEBUG:  {name} at_pos = {at_pos}")

            # Parse parameters
            params = self._parse_parameters(params_str)

            # # Parse rotation if present
            # rotation = [0, 0, 0]
            # rotation_match = re.search(r"ROTATED\s*\((.*?)\)\s*RELATIVE", block)
            # if rotation_match:
            #     rotation = [
            #         float(x.strip()) for x in rotation_match.group(1).split(",")
            #     ]
            # Parse rotation if present
            rotation = [0, 0, 0]
            rotation_match = re.search(r"ROTATED\s*\((.*?)\)\s*RELATIVE", block)
            if rotation_match:
                rotation = []
                for coord in rotation_match.group(1).split(","):
                    coord = coord.strip()
                    try:
                        # Try direct conversion first
                        rotation.append(float(coord))
                    except ValueError:
                        # Try to evaluate as expression with instrument parameters
                        try:
                            # # Create evaluation context with parameters and math functions
                            # import math

                            # eval_context = {
                            #     **self.parameters,
                            #     "PI": math.pi,
                            #     "sin": math.sin,
                            #     "cos": math.cos,
                            #     "tan": math.tan,
                            #     "atan": math.atan,
                            # }
                            value = eval(coord, {"__builtins__": {}}, eval_context)
                            rotation.append(float(value))
                            print(
                                f"DEBUG EVAL ROTATION: {name} - evaluated '{coord}' = {value}"
                            )
                        except Exception as e:
                            print(
                                f"WARNING: Could not evaluate rotation '{coord}' for {name}:  {e}, using 0.0"
                            )
                            rotation.append(0.0)

            component = McStasComponent(
                name, comp_type, at_pos, rel_to, rotation, params
            )
            self.components.append(component)

        # Calculate absolute positions
        self._calculate_absolute_positions()

        return self.components

    def _parse_parameters(self, params_str):
        """Parse component parameters"""
        params = {}
        if not params_str.strip():
            return params

        # Match parameter = value pairs
        param_matches = re.finditer(r"(\w+)\s*=\s*([^,\)]+)", params_str)
        for match in param_matches:
            key = match.group(1).strip()
            value_str = match.group(2).strip()
            # Try to convert to float, otherwise keep as string
            try:
                value = float(value_str)
            except ValueError:
                value = value_str.strip("\"'")
            params[key] = value

        return params

    # def _calculate_absolute_positions(self):
    #     """Calculate absolute positions for all components"""
    #     component_dict = {comp.name: comp for comp in self.components}

    #     def get_absolute_transform(comp):
    #         if comp.absolute_pos is not None:
    #             return comp.absolute_pos, comp.absolute_rotation

    #         if comp.relative_to == "ABSOLUTE":
    #             comp.absolute_pos = np.array(comp.at_pos)
    #             comp.absolute_rotation = np.array(comp.rotation)
    #         else:
    #             # Get parent's absolute position
    #             parent = component_dict.get(comp.relative_to)
    #             if parent:
    #                 parent_pos, parent_rot = get_absolute_transform(parent)
    #                 # For simplicity, just add positions (proper rotation transform would be more complex)
    #                 comp.absolute_pos = parent_pos + np.array(comp.at_pos)
    #                 comp.absolute_rotation = parent_rot + np.array(comp.rotation)
    #             else:
    #                 comp.absolute_pos = np.array(comp.at_pos)
    #                 comp.absolute_rotation = np.array(comp.rotation)

    #         return comp.absolute_pos, comp.absolute_rotation

    #     for comp in self.components:
    #         get_absolute_transform(comp)

    def _calculate_absolute_positions(self):
        """Calculate absolute positions for all components"""
        component_dict = {comp.name: comp for comp in self.components}

        def get_absolute_transform(comp):
            if comp.absolute_pos is not None:
                return comp.absolute_pos, comp.absolute_rotation

            if comp.relative_to == "ABSOLUTE":
                comp.absolute_pos = np.array(comp.at_pos)
                comp.absolute_rotation = np.array(comp.rotation)
                print(f"DEBUG POS: {comp.name} ABSOLUTE at {comp.absolute_pos}")
            else:
                # Get parent's absolute position
                parent = component_dict.get(comp.relative_to)
                if parent:
                    parent_pos, parent_rot = get_absolute_transform(parent)
                    # For simplicity, just add positions (proper rotation transform would be more complex)
                    comp.absolute_pos = parent_pos + np.array(comp.at_pos)
                    comp.absolute_rotation = parent_rot + np.array(comp.rotation)
                    print(
                        f"DEBUG POS: {comp.name} RELATIVE to {comp.relative_to}:  parent={parent_pos}, offset={comp.at_pos}, result={comp.absolute_pos}"
                    )
                else:
                    comp.absolute_pos = np.array(comp.at_pos)
                    comp.absolute_rotation = np.array(comp.rotation)
                    print(
                        f"DEBUG POS: {comp.name} no parent found for {comp.relative_to}, using {comp.absolute_pos}"
                    )

            # print(f"DEBUG ROT: {comp.name} rotation = {comp.absolute_rotation}")

            return comp.absolute_pos, comp.absolute_rotation

        for comp in self.components:
            get_absolute_transform(comp)


class McStasVisualizer:
    """Visualizes McStas instruments using pythreejs"""

    # Color scheme for different component types
    COMPONENT_COLORS = {
        "source": "#ff0000",  # Red
        "ess": "#ff0000",  # Red
        "guide": "#808080",  # Gray
        "slit": "#808080",  # Gray
        "monochromator": "#00ff00",  # Green
        "sample": "#ffff00",  # Yellow
        "sans": "#ffff00",  # Yellow
        "detector": "#0000ff",  # Blue
        "chopper": "#ff00ff",  # Magenta
        "slit": "#00ffff",  # Cyan
        "monitor": "#ffa500",  # Orange
        "psd": "#ffa500",  # Orange
        "arm": "#444444",  # Dark gray
        "cylinder": "#008000",
        "default": "#cccccc",  # Light gray
    }

    def __init__(self, components):
        self.components = components
        self.scene_objects = []

    def create_geometry(self, component):
        """Create pythreejs geometry for a component based on its type"""
        comp_type = component.component_type.lower()

        if comp_type in SKIP_TYPES:
            return None

        # print("component type", comp_type)
        params = component.params

        # Determine color based on component type
        color = self._get_color(comp_type)
        material = p3.MeshBasicMaterial(
            color=color,
            # transparent=True,
            # opacity=0.8,
            depthWrite=False,
            # side="DoubleSide",
            # depthTest=False
        )
        wireframe_material = p3.LineBasicMaterial(
            color="#000000", linewidth=1, opacity=0.8, transparent=True
        )

        geometry = None
        mesh = None

        print(f"Creating geometry for component: {component.name} ({comp_type})")

        # ESS_butterfly (source)
        if "ess" in comp_type or "source" in comp_type:
            yheight = params.get("yheight", 0.03)
            width = params.get("xwidth", 0.03)
            depth = 0.01
            geometry = p3.BoxGeometry(width, yheight, depth)
            mesh = p3.Mesh(geometry, material)

        # Guide, Slit
        elif "guide" in comp_type or "slit" in comp_type:
            width = params.get("xwidth", params.get("w1", 0.05))
            height = params.get("yheight", params.get("h1", 0.05))
            length = params.get("l", params.get("length", 0.5))
            geometry = p3.BoxGeometry(width, height, length)
            mesh = p3.Mesh(geometry, material)
        # Chopper
        elif "chopper" in comp_type or "diskchopper" in comp_type:
            # return None
            radius = params.get("radius", params.get("R", 0.3))
            thickness = params.get("thickness", 0.02)
            # geometry = p3.CylinderGeometry(radius, radius, thickness, 32)
            # mesh = p3.Mesh(geometry, material)

            # # Convert theta_0 to radians
            # theta_0 = params.get("theta_0", 0)
            # print("Chopper theta_0:", theta_0)
            # theta_rad = theta_0 * np.pi / 180

            # # Calculate how much of the disk to show (everything except the opening)
            # # The opening is centered around the beam path (at the calculated offset)
            # thetaLength = 2 * np.pi - theta_rad if theta_rad > 0 else 2 * np.pi
            # # Start angle offset so the opening is centered around the y-axis (beam path)
            # thetaStart = theta_rad / 2 if theta_rad > 0 else 0

            thetaStart = 0
            thetaLength = 2 * np.pi

            geometry = p3.CylinderGeometry(
                radius, radius, thickness, 32, 1, False, thetaStart, thetaLength
            )
            mesh = p3.Mesh(geometry, material)

            # Rotate to face beam direction (z-axis)
            # mesh.rotation = (np.pi / 2, 0, 0, "XYZ")
            mesh.rotateX(np.pi / 2)
            print("Chopper geometry created", mesh.rotation)

            # Offset the chopper so beam passes through the opening
            # Beam passes at y = 0.5 * (radius + radius - yheight) = 0.5 * (2*radius - yheight)
            yheight = params.get("yheight", -1)
            if yheight > 0:
                beam_offset = 0.5 * (2 * radius - yheight)
                # Move the chopper down so the beam passes through the opening
                mesh.position = (0, -beam_offset, 0)

        # Sample (sphere or box)
        elif "sample" in comp_type or "sans" in comp_type:
            print("Creating sample geometry")
            # if "radius" in params or "R" in params:
            #     radius = params.get("radius", params.get("R", 0.01))
            #     geometry = p3.SphereGeometry(radius, 16, 16)
            # else:
            width = params.get("xwidth", 0.02)
            height = params.get("yheight", 0.02)
            depth = params.get("zthick", params.get("zdepth", 0.002))
            geometry = p3.BoxGeometry(width, height, depth)
            # print("Sample geometry created", geometry)
            mesh = p3.Mesh(geometry, material)

        # Detector (cylinder or box)
        elif "detector" in comp_type or "monitor" in comp_type or "psd" in comp_type:
            if "radius" in params:
                radius = params.get("radius", 0.5)
                height = params.get("yheight", 1.0)
                geometry = p3.CylinderGeometry(radius, radius, height, 32)
                mesh = p3.Mesh(geometry, material)
            else:
                width = params.get("xwidth", 0.5)
                height = params.get("yheight", 1.0)
                depth = params.get("zdepth", 0.01)
                geometry = p3.BoxGeometry(width, height, depth)
                mesh = p3.Mesh(geometry, material)

        # Union components (cylinders and spheres)
        elif "union_cylinder" in comp_type:
            radius = params.get("radius", 0.03)
            height = params.get("yheight", 0.5)
            geometry = p3.CylinderGeometry(radius, radius, height, 32)
            mesh = p3.Mesh(geometry, material)

        elif "union_sphere" in comp_type:
            radius = params.get("radius", 0.01)
            geometry = p3.SphereGeometry(radius, 16, 16)
            mesh = p3.Mesh(geometry, material)

        elif "union_box" in comp_type:
            width = params.get("xwidth", 0.1)
            height = params.get("yheight", 0.1)
            depth = params.get("zdepth", 0.1)
            geometry = p3.BoxGeometry(width, height, depth)
            mesh = p3.Mesh(geometry, material)

        # Incoherent scatterer (box)
        elif "incoherent" in comp_type:
            width = params.get("xwidth", 0.05)
            height = params.get("yheight", 0.05)
            depth = params.get("zdepth", 0.008)
            geometry = p3.BoxGeometry(width, height, depth)
            mesh = p3.Mesh(geometry, material)

        # # Arm (small sphere as marker)
        # elif "arm" in comp_type or "init" in comp_type:
        #     geometry = p3.SphereGeometry(0.005, 8, 8)
        #     mesh = p3.Mesh(geometry, material)

        # Default:  small box
        else:
            geometry = p3.BoxGeometry(0.02, 0.02, 0.02)
            mesh = p3.Mesh(geometry, material)

        if mesh:
            # Add wireframe edges
            if geometry:
                edges = p3.EdgesGeometry(geometry)
                wireframe = p3.LineSegments(edges, wireframe_material, linewidth=2)
                if "chopper" in comp_type or "diskchopper" in comp_type:
                    # Rotate wireframe to match chopper orientation
                    wireframe.rotateX(np.pi / 2)
                    wireframe.position = mesh.position
                group = p3.Group()
                group.add(mesh)
                group.add(wireframe)
                return group
            return mesh

        return None

    def _get_color(self, comp_type):
        """Get color for component type"""
        for key, color in self.COMPONENT_COLORS.items():
            # print(key, color, comp_type, key in comp_type)
            if key in comp_type:
                return color
        return self.COMPONENT_COLORS["default"]

    def create_scene(self, width=900, height=600, show_axes=True, show_labels=True):
        """Create the complete 3D scene"""
        scene = p3.Scene(children=[])

        # # Add ambient lighting
        ambient = p3.AmbientLight(intensity=1.0)
        scene.add(ambient)

        # Add axes helper
        if show_axes:
            axes = p3.AxesHelper(size=1)
            scene.add(axes)

        # Add components
        for component in self.components:
            geom = self.create_geometry(component)
            if geom:
                # Set position
                pos = component.absolute_pos
                geom.position = tuple(pos)

                # # Set rotation (convert degrees to radians)
                # rot = component.absolute_rotation * np.pi / 180
                # # geom.rotation = (rot[0], rot[1], rot[2], "XYZ")
                # geom.rotateX(rot[0])
                # geom.rotateY(rot[1])
                # geom.rotateZ(rot[2])
                # # print(
                # #     "debuuuuuuuug",
                # #     component.name,
                # #     geom.rotation,
                # #     component.absolute_rotation,
                # # )

                scene.add(geom)
                # print(geom.rotation)
                self.scene_objects.append(geom)

                # # Add label (simplified - using small text would require TextGeometry)
                # if show_labels and component.component_type.lower() not in [
                #     "arm",
                #     "init",
                #     "stop",
                # ]:
                #     # Create a small marker sphere for labeled components
                #     label_geom = p3.SphereGeometry(0.01, 8, 8)
                #     label_material = p3.MeshBasicMaterial(color="#ffffff")
                #     label_mesh = p3.Mesh(label_geom, label_material)
                #     label_mesh.position = (pos[0], pos[1] + 0.05, pos[2])
                #     scene.add(label_mesh)

        # Add beam path (line connecting components along z-axis)
        self._add_beam_path(scene)

        # Create camera
        camera = p3.PerspectiveCamera(
            position=[5, 3, 10], aspect=width / height, fov=50, near=0.01, far=2000
        )
        camera.lookAt([0, 0, 2])

        # Create renderer with orbit controls
        controls = p3.OrbitControls(controlling=camera)
        renderer = p3.Renderer(
            camera=camera, scene=scene, controls=[controls], width=width, height=height
        )

        return renderer

    def _add_beam_path(self, scene):
        """Add a line showing the nominal beam path"""
        # Find components along the beam (sorted by z position)
        beam_components = [c for c in self.components if c.absolute_pos[2] >= 0]
        beam_components.sort(key=lambda c: c.absolute_pos[2])

        if len(beam_components) < 2:
            return

        # Create line geometry
        points = [
            tuple(c.absolute_pos)
            # for c in beam_components[: min(len(beam_components), 10)]
            for c in beam_components
        ]

        # Create dashed line
        line_material = p3.LineBasicMaterial(
            color="#ff0000", linewidth=2, opacity=0.3, transparent=True
        )
        line_geom = p3.BufferGeometry(
            attributes={"position": p3.BufferAttribute(array=points, normalized=False)}
        )
        line = p3.Line(geometry=line_geom, material=line_material)
        scene.add(line)

    # Add this method to the McStasVisualizer class:
    def create_component_navigator(self, renderer):
        """Create a dropdown widget to navigate to components"""

        # Create dropdown options (skip utility components)
        # skip_types = [
        #     "union_init",
        #     "union_stop",
        #     "union_make_material",
        #     "incoherent_process",
        #     "powder_process",
        # ]

        # self.component_mapping = {
        #     f"{comp.name} ({comp.component_type})": i
        #     for i, comp in enumerate(self.components)
        #     if comp.component_type.lower() not in skip_types
        # }

        # dropdown = ipw.Combobox(
        #     options=tuple(self.component_mapping.keys()),
        #     description="Take me to: ",
        #     style={"description_width": "60px"},
        # )

        component_options = [
            (f"{comp.name} ({comp.component_type})", i)
            for i, comp in enumerate(self.components)
            if comp.component_type.lower() not in SKIP_TYPES
        ]

        dropdown = ipw.Dropdown(
            options=component_options,
            description="Take me to: ",
            style={"description_width": "initial"},
        )

        # if not component_options:
        #     component_options = [
        #         (f"{comp.name} ({comp.component_type})", i)
        #         for i, comp in enumerate(self.components)
        #     ]

        def on_component_select(change):
            # idx = self.component_mapping.get(change["new"], None)
            # if idx is None:
            #     return
            idx = change["new"]
            component = self.components[idx]
            pos = component.absolute_pos

            # Calculate camera position (offset to top-left and back)
            # Determine a reasonable distance based on component type
            comp_type = component.component_type.lower()

            if (
                "source" in comp_type
                or "detector" in comp_type
                or "container" in comp_type
            ):
                distance = 2.0
            elif "sample" in comp_type:
                distance = 0.5
            else:
                distance = 1.0

            # Position camera at an angle (top-left-front)
            cam_x = pos[0] - distance * 0.5
            cam_y = pos[1] + distance * 0.7
            cam_z = pos[2] + distance * 0.5

            # Update camera position
            renderer.camera.position = [cam_x, cam_y, cam_z]

            # Update controls target to component center
            renderer.controls[0].target = [float(pos[0]), float(pos[1]), float(pos[2])]

            # Reset the camera to look at the target
            renderer.camera.lookAt([float(pos[0]), float(pos[1]), float(pos[2])])

        dropdown.observe(on_component_select, names="value")

        return dropdown

    # Modify the display method:
    def display(self, width=900, height=600, show_axes=True, show_labels=True):
        """Display the instrument in Jupyter with navigation controls"""
        renderer = self.create_scene(width, height, show_axes, show_labels)

        # Create component navigator dropdown
        navigator = self.create_component_navigator(renderer)

        # Display in a VBox
        # display(VBox([navigator, renderer]))

        return ipw.VBox([navigator, renderer])

    # def display(self, width=900, height=600, show_axes=True, show_labels=True):
    #     """Display the instrument in Jupyter"""
    #     renderer = self. create_scene(width, height, show_axes, show_labels)
    #     # display(renderer)
    #     return renderer


def visualize_mcstas_instrument(
    instrument_file, width=900, height=600, show_axes=True, show_labels=True
):
    """
    Main function to visualize a McStas instrument file

    Parameters:
    -----------
    instrument_file : str
        Path to the McStas instrument file
    width : int
        Width of the display in pixels
    height : int
        Height of the display in pixels
    show_axes : bool
        Whether to show coordinate axes
    show_labels : bool
        Whether to show component labels

    Returns:
    --------
    renderer : pythreejs.Renderer
        The renderer object that can be further manipulated
    """
    # Parse the instrument
    parser = McStasInstrumentParser(instrument_file)
    components = parser.parse()

    print(f"Parsed instrument:  {parser.instrument_name}")
    print(f"Found {len(components)} components")

    # Create and display visualization
    visualizer = McStasVisualizer(components)
    return visualizer.display(
        width=width, height=height, show_axes=show_axes, show_labels=show_labels
    )


# # Convenience function to list components
# def list_components(instrument_file):
#     """List all components in an instrument file"""
#     parser = McStasInstrumentParser(instrument_file)
#     components = parser.parse()

#     print(f"Instrument: {parser.instrument_name}")
#     print(f"\nComponents ({len(components)}):")
#     print("-" * 80)
#     print(f"{'Name':<20} {'Type':<25} {'Position (x,y,z)':<30}")
#     print("-" * 80)

#     for comp in components:
#         pos_str = f"({comp.absolute_pos[0]:.3f}, {comp.absolute_pos[1]:.3f}, {comp.absolute_pos[2]:.3f})"
#         print(f"{comp.name:<20} {comp.component_type:<25} {pos_str:<30}")

#     return components
