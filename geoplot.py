"""
geoplot.py - Geographic Visualization for Simulations
----------------------------------------------------

This module creates interactive 3D visualizations of simulation data over time.
It generates an HTML file using Cesium Ion for geographic rendering and creates
a GeoJSON file with the data needed for the visualization.

Example usage:
```python
from agent_torch.visualize import GeoPlot

# Create a visualizer with your configuration
viz = GeoPlot(config, {
    "cesium_token": "your-cesium-token",        # Cesium Ion API token
    "step_time": 3600,                          # Seconds between steps
    "coordinates": "agents/consumers/location",  # Path to coordinates in state
    "feature": "agents/consumers/money_spent",   # Path to values to visualize
    "visualization_type": "color",              # "color" or "size"
})

# In your simulation loop
for i in range(num_episodes):
    runner.step(steps_per_episode)
    viz.render(runner.state_trajectory)
```
"""

import re
import json
from typing import Dict, List, Any, Optional, Tuple

import pandas as pd
import numpy as np

from string import Template
from agent_torch.core.helpers import get_by_path

# HTML template for Cesium visualization (with improved readability and annotations)
geoplot_template = """
<!doctype html>
<html lang="en">
	<head>
		<meta charset="UTF-8" />
		<meta name="viewport" content="width=device-width, initial-scale=1.0" />
		<title>Geographic Visualization - $simulationName</title>
		<script src="https://cesium.com/downloads/cesiumjs/releases/1.95/Build/Cesium/Cesium.js"></script>
		<link href="https://cesium.com/downloads/cesiumjs/releases/1.95/Build/Cesium/Widgets/widgets.css" rel="stylesheet" />
		<style>
			body { margin: 0; padding: 0; overflow: hidden; }
			#cesiumContainer { width: 100%; height: 100vh; }
			.info-box { 
				position: absolute; top: 10px; left: 10px; 
				background: rgba(40, 40, 40, 0.7); color: white;
				padding: 10px; border-radius: 4px; z-index: 1000;
				font-family: Arial, sans-serif; max-width: 250px;
			}
		</style>
	</head>
	<body>
		<div id="cesiumContainer"></div>
		<div class="info-box">
			<h3>$simulationName</h3>
			<p>Visualizing: <strong>$featureName</strong></p>
			<p>Using: <strong>$visualType</strong> mapping</p>
		</div>
		<script>
			// Initialize Cesium with your token
			Cesium.Ion.defaultAccessToken = '$accessToken';
			const viewer = new Cesium.Viewer('cesiumContainer');
			
			// Color interpolation function
			function mixColors(color1, color2, factor) {
				const result = new Cesium.Color();
				result.red = color1.red + factor * (color2.red - color1.red);
				result.green = color1.green + factor * (color2.green - color1.green);
				result.blue = color1.blue + factor * (color2.blue - color1.blue);
				result.alpha = '$visualType' === 'size' ? 0.2 : 
					color1.alpha + factor * (color2.alpha - color1.alpha);
				return result;
			}
			
			// Map value to color (blue to red)
			function getValueColor(value, min, max) {
				const factor = (value - min) / (max - min);
				return mixColors(Cesium.Color.BLUE, Cesium.Color.RED, factor);
			}
			
			// Map value to point size
			function getValueSize(value, min, max) {
				const factor = (value - min) / (max - min);
				return 100 * (1 + factor);
			}
			
			// Process GeoJSON data into time series format
			function processTimeSeriesData(geoJson) {
				const timeSeriesMap = new Map();
				let minValue = Infinity;
				let maxValue = -Infinity;

				geoJson.features.forEach((feature) => {
					const id = feature.properties.id;
					const time = Cesium.JulianDate.fromIso8601(feature.properties.time);
					const value = feature.properties.value;
					const coordinates = feature.geometry.coordinates;

					if (!timeSeriesMap.has(id)) {
						timeSeriesMap.set(id, []);
					}
					timeSeriesMap.get(id).push({ time, value, coordinates });

					minValue = Math.min(minValue, value);
					maxValue = Math.max(maxValue, value);
				});

				return { timeSeriesMap, minValue, maxValue };
			}
			
			// Create Cesium entities for visualization
			function createEntities(timeSeriesData, startTime, endTime) {
				const dataSource = new Cesium.CustomDataSource('Simulation Data');
				
				for (const [id, timeSeries] of timeSeriesData.timeSeriesMap) {
					const entity = new Cesium.Entity({
						id: id,
						name: `Entity ${id}`,
						availability: new Cesium.TimeIntervalCollection([
							new Cesium.TimeInterval({
								start: startTime,
								stop: endTime,
							}),
						]),
						position: new Cesium.SampledPositionProperty(),
						point: {
							pixelSize: '$visualType' === 'size' ? 
								new Cesium.SampledProperty(Number) : 10,
							color: new Cesium.SampledProperty(Cesium.Color),
							outlineColor: Cesium.Color.BLACK,
							outlineWidth: 1,
						},
						properties: {
							value: new Cesium.SampledProperty(Number),
						},
					});

					// Add data samples for each time point
					timeSeries.forEach(({ time, value, coordinates }) => {
						const position = Cesium.Cartesian3.fromDegrees(
							coordinates[0],  // longitude
							coordinates[1]   // latitude
						);
						
						entity.position.addSample(time, position);
						entity.properties.value.addSample(time, value);
						entity.point.color.addSample(
							time,
							getValueColor(value, timeSeriesData.minValue, timeSeriesData.maxValue)
						);

						if ('$visualType' === 'size') {
							entity.point.pixelSize.addSample(
								time,
								getValueSize(value, timeSeriesData.minValue, timeSeriesData.maxValue)
							);
						}
					});

					dataSource.entities.add(entity);
				}
				return dataSource;
			}

			// Load the data and set up the visualization
			const datasets = $data;
			const start = Cesium.JulianDate.fromIso8601('$startTime');
			const stop = Cesium.JulianDate.fromIso8601('$stopTime');
			
			// Configure clock settings
			viewer.clock.startTime = start.clone();
			viewer.clock.stopTime = stop.clone();
			viewer.clock.currentTime = start.clone();
			viewer.clock.clockRange = Cesium.ClockRange.LOOP_STOP;
			viewer.clock.multiplier = 3600;  // 1 hour per second
			
			viewer.timeline.zoomTo(start, stop);
			
			// Process each dataset
			for (const geoJsonData of datasets) {
				const timeSeriesData = processTimeSeriesData(geoJsonData);
				const dataSource = createEntities(timeSeriesData, start, stop);
				viewer.dataSources.add(dataSource);
				viewer.zoomTo(dataSource);
			}
		</script>
	</body>
</html>
"""


def get_data_from_path(state: Dict, path: str) -> Any:
    """
    Retrieves data from a nested state dictionary using a path string.
    
    Args:
        state: The state dictionary with nested data
        path: A string path with '/' separators (e.g. 'agents/consumers/money')
    
    Returns:
        The value found at the specified path
    """
    path_parts = re.split("/", path)
    return get_by_path(state, path_parts)


class GeoPlot:
    """
    Creates interactive geographic visualizations of simulation data.
    
    This class takes simulation state data and produces:
    1. A GeoJSON file with all data points
    2. An HTML file with a Cesium-based 3D visualization
    """
    
    def __init__(self, config: Dict[str, Any], options: Dict[str, Any]):
        """
        Set up the visualization tool with configuration settings.
        
        Args:
            config: Simulation configuration dictionary
            options: Visualization options with these keys:
                - cesium_token: API token for Cesium Ion
                - step_time: Seconds between simulation steps
                - coordinates: Path to coordinates in state
                - feature: Path to values to visualize
                - visualization_type: "color" or "size" (optional)
        """
        self.config = config
        
        # Extract required options
        self.cesium_token = options["cesium_token"]
        self.step_time = options["step_time"]
        self.coordinates_path = options["coordinates"]
        self.feature_path = options["feature"]
        
        # Optional settings with defaults
        self.visualization_type = options.get("visualization_type", "color")
        
        # Validate visualization type
        if self.visualization_type not in ["color", "size"]:
            print(f"Warning: Unknown visualization type '{self.visualization_type}', using 'color'")
            self.visualization_type = "color"
        
        # Extract feature name for display
        self.feature_name = self.feature_path.split("/")[-1].replace("_", " ").title()
        
        # Simulation settings
        self.simulation_name = config["simulation_metadata"]["name"]
        self.num_episodes = config["simulation_metadata"]["num_episodes"]
        self.steps_per_episode = config["simulation_metadata"]["num_steps_per_episode"]
        
        print(f"GeoPlot initialized to visualize '{self.feature_name}' using {self.visualization_type}")

    def render(self, state_trajectory: List) -> Tuple[str, str]:
        """
        Generate visualization files from simulation data.
        
        Args:
            state_trajectory: List of simulation states over time
        
        Returns:
            Tuple of (geojson_path, html_path) with the output file paths
        """
        # Initialize empty data containers
        coordinates = []
        values = []
        
        # Set output file paths
        geojson_path = f"{self.simulation_name}.geojson"
        html_path = f"{self.simulation_name}.html"
        
        print(f"Creating visualization for '{self.simulation_name}'...")
        
        # Extract data from each state
        for i, step_states in enumerate(state_trajectory[:-1]):
            # Get final state of the step
            final_state = step_states[-1]
            
            try:
                # Get coordinates and values, convert to lists if needed
                coords = get_data_from_path(final_state, self.coordinates_path)
                if isinstance(coords, np.ndarray):
                    coords = coords.tolist()
                coordinates = coords  # Store latest coordinates
                
                # Get values for this step
                vals = get_data_from_path(final_state, self.feature_path)
                if isinstance(vals, np.ndarray):
                    vals = vals.flatten().tolist()
                values.append(vals)
                
            except Exception as e:
                print(f"Error extracting data (step {i}): {e}")
                print(f"Check paths: coordinates='{self.coordinates_path}', feature='{self.feature_path}'")
                return None, None
        
        # Generate timestamps for each step
        start_time = pd.Timestamp.utcnow()
        total_steps = self.num_episodes * self.steps_per_episode
        
        timestamps = [
            start_time + pd.Timedelta(seconds=i * self.step_time)
            for i in range(total_steps)
        ]
        
        # Create GeoJSON collection
        geojson_collection = self._create_geojson(coordinates, values, timestamps)
        
        # Write GeoJSON to file
        try:
            with open(geojson_path, "w", encoding="utf-8") as f:
                json.dump(geojson_collection, f, ensure_ascii=False, indent=2)
            print(f"GeoJSON data written to {geojson_path}")
        except Exception as e:
            print(f"Error writing GeoJSON: {e}")
            return None, None
        
        # Create HTML visualization
        try:
            html_content = self._create_html(geojson_collection, timestamps)
            with open(html_path, "w", encoding="utf-8") as f:
                f.write(html_content)
            print(f"Visualization created at {html_path}")
        except Exception as e:
            print(f"Error creating HTML: {e}")
            return None, None
            
        return geojson_path, html_path
    
    def _create_geojson(self, coordinates: List, values: List, timestamps: List) -> List:
        """
        Create GeoJSON data structure from coordinates and values.
        
        Args:
            coordinates: List of coordinate pairs [[lat, lon], ...]
            values: List of value lists for each time step
            timestamps: List of timestamps for each step
            
        Returns:
            List of GeoJSON FeatureCollection objects
        """
        geojson_collection = []
        
        # For each entity (coordinate pair)
        for entity_idx, coord in enumerate(coordinates):
            features = []
            
            # For each time step with values
            for time_idx, (time, value_list) in enumerate(zip(timestamps, values)):
                # Make sure we have a value for this entity
                if entity_idx < len(value_list):
                    # Create a feature for this entity at this time
                    # Note: GeoJSON uses [longitude, latitude] order
                    features.append({
                        "type": "Feature",
                        "geometry": {
                            "type": "Point",
                            "coordinates": [coord[1], coord[0]]  # [lon, lat]
                        },
                        "properties": {
                            "id": f"entity_{entity_idx}",
                            "value": value_list[entity_idx],
                            "time": time.isoformat()
                        }
                    })
            
            # Add this entity's features if we have any
            if features:
                geojson_collection.append({
                    "type": "FeatureCollection",
                    "features": features
                })
                
        return geojson_collection
    
    def _create_html(self, geojson_data: List, timestamps: List) -> str:
        """
        Create HTML visualization from template and data.
        
        Args:
            geojson_data: GeoJSON data structure
            timestamps: List of timestamps for visualization
            
        Returns:
            HTML content as string
        """
        # Create visualization from template
        template = Template(geoplot_template)
        
        # Replace template variables with actual values
        return template.substitute({
            "simulationName": self.simulation_name,
            "featureName": self.feature_name,
            "accessToken": self.cesium_token,
            "startTime": timestamps[0].isoformat(),
            "stopTime": timestamps[-1].isoformat(),
            "data": json.dumps(geojson_data),
            "visualType": self.visualization_type
        })