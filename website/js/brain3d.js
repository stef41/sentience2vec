/**
 * Brain3D Module - BrainBrowser Surface Viewer Integration
 * Visualizes receptor density, drug binding, and effects on 3D brain
 */

(function() {
    'use strict';

    // Brain region to approximate 3D coordinates mapping
    // Based on MNI coordinates for fsaverage brain
    var REGION_COORDS = {
        // Frontal regions
        'prefrontal_cortex': { x: 0, y: 60, z: 20 },
        'prefrontal_cortex_layer_ii': { x: 0, y: 60, z: 20 },
        'prefrontal_cortex_layer_v': { x: 0, y: 55, z: 15 },
        'orbitofrontal_cortex': { x: 0, y: 50, z: -15 },
        'anterior_cingulate_cortex': { x: 0, y: 35, z: 25 },
        'frontal_cortex': { x: 0, y: 45, z: 30 },
        
        // Motor regions
        'motor_cortex': { x: 0, y: 0, z: 60 },
        'premotor_cortex': { x: 0, y: 15, z: 55 },
        'supplementary_motor_area': { x: 0, y: 5, z: 60 },
        
        // Temporal regions
        'temporal_cortex': { x: 55, y: -15, z: -10 },
        'superior_temporal_gyrus': { x: 55, y: -20, z: 5 },
        'inferior_temporal_cortex': { x: 50, y: -30, z: -20 },
        'auditory_cortex': { x: 55, y: -25, z: 10 },
        
        // Parietal regions
        'parietal_cortex': { x: 0, y: -45, z: 50 },
        'somatosensory_cortex': { x: 0, y: -25, z: 55 },
        'posterior_parietal_cortex': { x: 0, y: -55, z: 45 },
        
        // Occipital regions
        'occipital_cortex': { x: 0, y: -85, z: 10 },
        'visual_cortex': { x: 0, y: -90, z: 5 },
        'visual_cortex_layer_iv': { x: 0, y: -90, z: 5 },
        
        // Limbic regions
        'hippocampus': { x: 25, y: -25, z: -15 },
        'hippocampus_ca1': { x: 25, y: -25, z: -15 },
        'hippocampus_ca3': { x: 28, y: -22, z: -15 },
        'hippocampus_dentate_gyrus': { x: 22, y: -28, z: -18 },
        'amygdala': { x: 22, y: -5, z: -20 },
        'amygdala_basolateral': { x: 22, y: -5, z: -20 },
        'amygdala_central': { x: 20, y: -3, z: -18 },
        'entorhinal_cortex': { x: 25, y: -10, z: -30 },
        'entorhinal_cortex_layer_ii': { x: 25, y: -10, z: -30 },
        'cingulate_cortex': { x: 0, y: 0, z: 35 },
        'posterior_cingulate': { x: 0, y: -45, z: 30 },
        'insula': { x: 38, y: 5, z: 0 },
        'insular_cortex': { x: 38, y: 5, z: 0 },
        
        // Subcortical
        'thalamus': { x: 0, y: -15, z: 5 },
        'thalamus_mediodorsal': { x: 5, y: -15, z: 8 },
        'hypothalamus': { x: 0, y: -5, z: -10 },
        'hypothalamus_paraventricular': { x: 3, y: -3, z: -8 },
        'striatum': { x: 15, y: 10, z: 5 },
        'caudate': { x: 12, y: 15, z: 10 },
        'putamen': { x: 25, y: 5, z: 5 },
        'nucleus_accumbens': { x: 10, y: 12, z: -8 },
        'globus_pallidus': { x: 18, y: 0, z: 0 },
        'ventral_tegmental_area': { x: 0, y: -20, z: -15 },
        'substantia_nigra': { x: 10, y: -20, z: -12 },
        'substantia_nigra_compacta': { x: 10, y: -20, z: -12 },
        
        // Brainstem
        'raphe_nuclei': { x: 0, y: -30, z: -25 },
        'dorsal_raphe_nucleus': { x: 0, y: -28, z: -22 },
        'median_raphe_nucleus': { x: 0, y: -32, z: -28 },
        'locus_coeruleus': { x: 5, y: -35, z: -25 },
        'periaqueductal_gray': { x: 0, y: -30, z: -20 },
        'brainstem': { x: 0, y: -35, z: -30 },
        
        // Cerebellum
        'cerebellum': { x: 0, y: -65, z: -35 },
        'cerebellar_cortex': { x: 25, y: -70, z: -30 },
        
        // Other
        'claustrum': { x: 32, y: 5, z: 5 },
        'septum': { x: 0, y: 8, z: 0 },
        'septum_lateral': { x: 5, y: 8, z: 0 },
        'basal_forebrain': { x: 0, y: 5, z: -10 },
        'olfactory_bulb': { x: 0, y: 35, z: -20 },
        'spinal_cord': { x: 0, y: -45, z: -50 }
    };

    // IUPHAR/BPS Guide to Pharmacology receptor IDs for source linking
    var RECEPTOR_GUIDE_IDS = {
        '5-HT1A': 1,
        '5-HT1B': 2,
        '5-HT1D': 3,
        '5-HT1E': 4,
        '5-HT1F': 5,
        '5-HT2A': 6,
        '5-HT2B': 7,
        '5-HT2C': 8,
        '5-HT3': 9,
        '5-HT4': 10,
        '5-HT5A': 11,
        '5-HT6': 12,
        '5-HT7': 13,
        'D1': 214,
        'D2': 215,
        'D3': 216,
        'D4': 217,
        'D5': 218,
        'MOR': 319,
        'DOR': 318,
        'KOR': 320,
        'NOP': 321,
        'CB1': 56,
        'CB2': 57,
        'NMDA': 459,
        'GABA-A': 408,
        'GABA_A_alpha1': 408,
        'GABA_A_alpha2': 408,
        'M1': 13,
        'M2': 14,
        'M3': 15,
        'M4': 16,
        'M5': 17,
        'alpha1A': 22,
        'alpha2A': 25,
        'alpha2B': 26,
        'alpha2C': 27,
        'beta1': 28,
        'beta2': 29,
        'H1': 262,
        'H2': 263,
        'H3': 264,
        'H4': 265,
        'SERT': 239,
        'DAT': 238,
        'NET': 240,
        'VMAT2': 245,
        'Sigma1': 377,
        'TAAR1': 365
    };

    function getGuideToPharmId(receptorName) {
        return RECEPTOR_GUIDE_IDS[receptorName] || '';
    }

    // Configure BrainBrowser
    BrainBrowser.config.set("worker_dir", "js/brainbrowser/workers/");

    // Store original colors for reset
    var originalColors = null;
    var currentReceptor = null;
    var currentEffect = null;
    var regionMarkers = []; // Store 3D markers for brain regions
    var showLabels = true;

    // Key brain regions to always display (major areas)
    var KEY_REGIONS = [
        'prefrontal_cortex',
        'hippocampus',
        'amygdala',
        'thalamus',
        'striatum',
        'cerebellum',
        'visual_cortex',
        'motor_cortex',
        'insula',
        'nucleus_accumbens',
        'raphe_nuclei',
        'locus_coeruleus'
    ];

    // Wait for DOM ready
    $(document).ready(function() {
        initBrainViewer();
        // Also populate selects immediately if DATA is available
        populateEffectSelect();
        populateReceptorSelect();
        populateDrugSelect();
    });

    function initBrainViewer() {
        var loadingDiv = document.getElementById('loading-indicator');
        
        if (!BrainBrowser.WEBGL_ENABLED) {
            document.getElementById('brainbrowser').innerHTML = 
                '<div class="error-message">WebGL is not supported. Please use Chrome, Firefox, or Edge.</div>';
            if (loadingDiv) loadingDiv.style.display = 'none';
            return;
        }

        window.brainViewer = BrainBrowser.SurfaceViewer.start("brainbrowser", function(viewer) {
            window.viewer = viewer;
            
            viewer.addEventListener("error", function(error) {
                console.error("BrainBrowser error:", error);
                if (loadingDiv) {
                    loadingDiv.textContent = "Error loading brain model";
                    loadingDiv.style.color = "#ff6b6b";
                }
            });
            
            viewer.addEventListener("displaymodel", function(event) {
                console.log("Model displayed:", event);
                if (loadingDiv) loadingDiv.style.display = 'none';
                
                viewer.setView("lateral");
                viewer.setClearColor(0x0c0c0c);
                
                // Store original vertex colors
                storeOriginalColors(viewer);
                
                // Skip 3D markers - we now use parcellation colors on the surface
                // createRegionMarkers(viewer);
                
                // Setup hover effects
                setupHoverEffects(viewer);
                
                // Populate UI
                populateEffectSelect();
                populateReceptorSelect();
                populateDrugSelect();
            });

            viewer.render();
            viewer.setClearColor(0x0c0c0c);
            loadBrainModel(viewer);
            setupControls(viewer);
            
            // Handle window resize
            window.addEventListener('resize', function() {
                if (viewer && viewer.updateViewport) {
                    viewer.updateViewport();
                }
            });
        });
    }

    // Store parcellation data
    var parcellationData = null;

    function loadBrainModel(viewer) {
        var loadingDiv = document.getElementById('loading-indicator');
        if (loadingDiv) loadingDiv.textContent = 'Loading brain surface...';
        
        // First load the parcellated brain JSON to get region data, then load model
        fetch("models/brain_parcellated.json")
            .then(function(response) { return response.json(); })
            .then(function(data) {
                parcellationData = data.parcellation;
                console.log("Parcellation loaded:", parcellationData.regions.length, "regions");
                
                // Now load the model after parcellation is ready
                viewer.loadModelFromURL("models/brain_parcellated.json", {
                    format: "json",
                    complete: function() {
                        console.log("Brain model loaded successfully");
                        if (loadingDiv) loadingDiv.style.display = 'none';
                        
                        // Apply default parcellation colors
                        if (parcellationData) {
                            applyParcellationColors(viewer);
                        }
                    },
                    error: function(error) {
                        console.error("Failed to load brain model:", error);
                        if (loadingDiv) {
                            loadingDiv.textContent = "Failed to load brain model";
                            loadingDiv.style.color = "#ff6b6b";
                        }
                    }
                });
            })
            .catch(function(err) {
                console.warn("Could not load parcellation:", err);
                // Still try to load model even if parcellation fails
                viewer.loadModelFromURL("models/brain_parcellated.json", {
                    format: "json",
                    complete: function() {
                        console.log("Brain model loaded (no parcellation)");
                        if (loadingDiv) loadingDiv.style.display = 'none';
                    },
                    error: function(error) {
                        console.error("Failed to load brain model:", error);
                        if (loadingDiv) {
                            loadingDiv.textContent = "Failed to load brain model";
                            loadingDiv.style.color = "#ff6b6b";
                        }
                    }
                });
            });
    }

    function storeOriginalColors(viewer) {
        var model = viewer.model;
        if (!model || !model.children || model.children.length === 0) return;
        
        var shape = model.children[0];
        if (shape && shape.geometry && shape.geometry.attributes && shape.geometry.attributes.color) {
            originalColors = new Float32Array(shape.geometry.attributes.color.array);
        }
    }

    // Hover state
    var hoveredRegion = null;
    var hoveredVertices = [];
    var raycaster = null;
    var mouse = null;

    function setupHoverEffects(viewer) {
        var container = document.getElementById('brainbrowser');
        if (!container) return;
        
        // Create tooltip element
        var tooltip = document.createElement('div');
        tooltip.id = 'brain-tooltip';
        tooltip.className = 'brain-tooltip';
        container.appendChild(tooltip);
        
        // Store mouse position for picking
        var mouseX = 0;
        var mouseY = 0;
        
        // Mouse move handler
        container.addEventListener('mousemove', function(event) {
            var rect = container.getBoundingClientRect();
            mouseX = event.clientX - rect.left;
            mouseY = event.clientY - rect.top;
            
            // Update tooltip position
            tooltip.style.left = (mouseX + 15) + 'px';
            tooltip.style.top = (mouseY + 15) + 'px';
            
            // Use BrainBrowser's built-in pick function
            handleHover(viewer, tooltip, mouseX, mouseY);
        });
        
        // Mouse leave handler
        container.addEventListener('mouseleave', function() {
            tooltip.style.display = 'none';
            clearHoverHighlight(viewer);
        });
        
        // Click handler for selecting region
        container.addEventListener('click', function(event) {
            if (hoveredRegion !== null && parcellationData) {
                showRegionDetails(hoveredRegion);
            }
        });
    }

    function handleHover(viewer, tooltip, mouseX, mouseY) {
        // Use BrainBrowser's built-in pick function
        var pickResult = viewer.pick(mouseX, mouseY);
        
        if (pickResult && pickResult.index !== undefined && parcellationData && parcellationData.vertex_labels && parcellationData.regions) {
            var vertexIndex = pickResult.index;
            var labels = parcellationData.vertex_labels;
            
            if (vertexIndex < labels.length) {
                var regionIndex = labels[vertexIndex];
                var regionName = parcellationData.regions[regionIndex] || 'Unknown';
                
                if (regionIndex !== hoveredRegion) {
                    // Clear previous highlight
                    clearHoverHighlight(viewer);
                    
                    // Highlight new region
                    hoveredRegion = regionIndex;
                    highlightRegion(viewer, regionIndex);
                }
                
                // Update tooltip
                var displayName = formatRegionName(regionName);
                var tooltipContent = '<strong>' + displayName + '</strong>';
                
                // Add receptor info if we're in receptor view mode
                if (currentReceptor && window.DATA) {
                    var receptor = window.DATA.receptors[currentReceptor];
                    if (receptor && receptor.regions) {
                        // Use the atlas index directly to find matching data regions
                        var regionData = findRegionDensityByAtlasIndex(receptor, regionIndex);
                        if (regionData) {
                            tooltipContent += '<br><span class="tooltip-receptor">' + currentReceptor + '</span>';
                            tooltipContent += '<br><span class="tooltip-density">Density: ' + regionData.density + ' fmol/mg</span>';
                            tooltipContent += '<br><span class="tooltip-region-match">(' + regionData.region.replace(/_/g, ' ') + ')</span>';
                        } else {
                            tooltipContent += '<br><span class="tooltip-receptor">' + currentReceptor + '</span>';
                            tooltipContent += '<br><span class="tooltip-no-data">No density data</span>';
                        }
                    }
                }
                
                // Add effect info if we're in effect view mode
                if (currentEffect && window.DATA && window.DATA.effects) {
                    var effect = window.DATA.effects[currentEffect];
                    if (effect) {
                        var effectDisplay = effect.name || currentEffect.replace(/_/g, ' ');
                        tooltipContent += '<br><span class="tooltip-effect">Effect: ' + effectDisplay + '</span>';
                    }
                }
                
                tooltip.innerHTML = tooltipContent;
                tooltip.style.display = 'block';
            } else {
                tooltip.style.display = 'none';
            }
        } else {
            tooltip.style.display = 'none';
            if (hoveredRegion !== null) {
                clearHoverHighlight(viewer);
                hoveredRegion = null;
            }
        }
        
        viewer.updated = true;
    }

    function highlightRegion(viewer, regionIndex) {
        if (!parcellationData) return;
        
        var model = viewer.model;
        if (!model || !model.children || model.children.length === 0) return;
        
        var shape = model.children[0];
        if (!shape || !shape.geometry) return;
        
        var colors = shape.geometry.attributes.color;
        if (!colors) return;
        
        var colorArray = colors.array;
        var labels = parcellationData.vertex_labels;
        
        hoveredVertices = [];
        
        // Brighten vertices in this region
        for (var i = 0; i < labels.length; i++) {
            if (labels[i] === regionIndex) {
                hoveredVertices.push({
                    index: i,
                    r: colorArray[i * 4],
                    g: colorArray[i * 4 + 1],
                    b: colorArray[i * 4 + 2]
                });
                
                // Brighten - add glow effect
                colorArray[i * 4] = Math.min(1, colorArray[i * 4] + 0.4);
                colorArray[i * 4 + 1] = Math.min(1, colorArray[i * 4 + 1] + 0.4);
                colorArray[i * 4 + 2] = Math.min(1, colorArray[i * 4 + 2] + 0.4);
            }
        }
        
        colors.needsUpdate = true;
    }

    function clearHoverHighlight(viewer) {
        if (hoveredVertices.length === 0) return;
        
        var model = viewer.model;
        if (!model || !model.children || model.children.length === 0) return;
        
        var shape = model.children[0];
        if (!shape || !shape.geometry) return;
        
        var colors = shape.geometry.attributes.color;
        if (!colors) return;
        
        var colorArray = colors.array;
        
        // Restore original colors
        hoveredVertices.forEach(function(v) {
            colorArray[v.index * 4] = v.r;
            colorArray[v.index * 4 + 1] = v.g;
            colorArray[v.index * 4 + 2] = v.b;
        });
        
        hoveredVertices = [];
        colors.needsUpdate = true;
        hoveredRegion = null;
    }

    function formatRegionName(name) {
        // Convert Destrieux atlas names to readable format
        return name
            .replace(/^G_/, 'Gyrus: ')
            .replace(/^S_/, 'Sulcus: ')
            .replace(/_/g, ' ')
            .replace(/-/g, ' ')
            .replace(/\b\w/g, function(l) { return l.toUpperCase(); });
    }

    function findMatchingKGRegion(atlasRegionName) {
        // This function is no longer used for tooltip density lookup.
        // Use getDataRegionsForAtlasIndex instead.
        return null;
    }
    
    function findRegionDensity(receptor, kgRegionPrefixes) {
        // This is the old method - kept for backwards compatibility but
        // the new method uses atlas index directly
        if (!receptor || !receptor.regions || !kgRegionPrefixes) return null;
        
        for (var i = 0; i < kgRegionPrefixes.length; i++) {
            var prefix = kgRegionPrefixes[i];
            var found = receptor.regions.find(function(r) {
                return r.region.toLowerCase().includes(prefix);
            });
            if (found) return found;
        }
        return null;
    }
    
    // Build reverse mapping: atlas index -> list of data region names
    function getDataRegionsForAtlasIndex(atlasIndex) {
        var dataRegions = [];
        Object.keys(REGION_TO_ATLAS).forEach(function(dataRegion) {
            var atlasIndices = REGION_TO_ATLAS[dataRegion];
            if (atlasIndices && atlasIndices.indexOf(atlasIndex) !== -1) {
                dataRegions.push(dataRegion);
            }
        });
        return dataRegions;
    }
    
    function findRegionDensityByAtlasIndex(receptor, atlasIndex) {
        if (!receptor || !receptor.regions) return null;
        
        var dataRegions = getDataRegionsForAtlasIndex(atlasIndex);
        
        // Find the highest density among all matching data regions
        var bestMatch = null;
        receptor.regions.forEach(function(regionData) {
            if (dataRegions.indexOf(regionData.region) !== -1) {
                if (!bestMatch || regionData.density > bestMatch.density) {
                    bestMatch = regionData;
                }
            }
        });
        
        return bestMatch;
    }

    function showRegionDetails(regionIndex) {
        if (!parcellationData) return;
        
        var regionName = parcellationData.regions[regionIndex];
        var displayName = formatRegionName(regionName);
        
        // Find receptors that have density in this region using atlas index
        var receptorsInRegion = [];
        if (window.DATA) {
            Object.keys(window.DATA.receptors).forEach(function(receptorName) {
                var receptor = window.DATA.receptors[receptorName];
                var found = findRegionDensityByAtlasIndex(receptor, regionIndex);
                if (found) {
                    receptorsInRegion.push({
                        name: receptorName,
                        density: found.density
                    });
                }
            });
        }
        
        // Sort by density
        receptorsInRegion.sort(function(a, b) { return b.density - a.density; });
        
        // Update info panel
        var regionInfo = document.getElementById('region-info');
        if (regionInfo) {
            var html = '<div class="receptor-info">';
            html += '<h4>📍 ' + displayName + '</h4>';
            
            if (receptorsInRegion.length > 0) {
                html += '<div class="density-regions">';
                html += '<h5>🧬 Receptors in this region</h5>';
                html += '<ul class="region-list">';
                
                var maxDensity = receptorsInRegion[0].density;
                receptorsInRegion.slice(0, 8).forEach(function(r) {
                    var barWidth = (r.density / maxDensity * 100);
                    html += '<li>';
                    html += '<span class="region-name">' + r.name + '</span>';
                    html += '<div class="density-bar-container">';
                    html += '<div class="density-bar" style="width:' + barWidth + '%"></div>';
                    html += '</div>';
                    html += '<span class="density-value">' + r.density + '</span>';
                    html += '</li>';
                });
                
                html += '</ul></div>';
            } else {
                html += '<p class="placeholder">No receptor data for this region</p>';
            }
            
            html += '</div>';
            regionInfo.innerHTML = html;
        }
    }

    function createRegionMarkers(viewer) {
        var THREE = BrainBrowser.SurfaceViewer.THREE;
        if (!THREE) {
            console.error("THREE not available");
            return;
        }
        
        // The scene in BrainBrowser is viewer.model.scene or we add to model itself
        var scene = viewer.model;
        if (!scene) {
            console.error("viewer.model not available");
            return;
        }
        
        // Remove any existing markers
        regionMarkers.forEach(function(marker) {
            if (scene.remove) scene.remove(marker);
        });
        regionMarkers = [];
        
        // Create markers for key regions
        KEY_REGIONS.forEach(function(regionName) {
            var coords = REGION_COORDS[regionName];
            if (!coords) return;
            
            try {
                // Create a group for marker + label
                var markerGroup = new THREE.Object3D();
                markerGroup.position.set(coords.x, coords.y, coords.z);
                markerGroup.userData.regionName = regionName;
                
                // Create small glowing sphere
                var sphereGeo = new THREE.SphereGeometry(2, 16, 16);
                var sphereMat = new THREE.MeshBasicMaterial({ 
                    color: 0x00ffaa, 
                    transparent: true, 
                    opacity: 0.8 
                });
                var sphere = new THREE.Mesh(sphereGeo, sphereMat);
                markerGroup.add(sphere);
                
                // Create outer glow ring - use number for DoubleSide if constant not available
                var ringGeo = new THREE.RingGeometry(2.5, 4, 32);
                var ringMat = new THREE.MeshBasicMaterial({ 
                    color: 0x00ffaa, 
                    transparent: true, 
                    opacity: 0.3,
                    side: THREE.DoubleSide || 2
                });
                var ring = new THREE.Mesh(ringGeo, ringMat);
                if (viewer.camera) ring.lookAt(viewer.camera.position);
                markerGroup.add(ring);
                
                scene.add(markerGroup);
                regionMarkers.push(markerGroup);
            } catch (e) {
                console.error("Error creating marker for", regionName, e);
            }
        });
        
        // Create HTML overlay for labels
        createLabelOverlay();
        viewer.updated = true;
        console.log("Created", regionMarkers.length, "region markers");
    }

    function createLabelOverlay() {
        // Remove existing overlay
        var existingOverlay = document.getElementById('region-labels-overlay');
        if (existingOverlay) existingOverlay.remove();
        
        var overlay = document.createElement('div');
        overlay.id = 'region-labels-overlay';
        overlay.style.cssText = 'position:absolute;top:0;left:0;width:100%;height:100%;pointer-events:none;overflow:hidden;';
        
        KEY_REGIONS.forEach(function(regionName) {
            var label = document.createElement('div');
            label.className = 'region-label';
            label.id = 'label-' + regionName;
            label.textContent = formatRegionName(regionName);
            overlay.appendChild(label);
        });
        
        document.getElementById('brainbrowser').appendChild(overlay);
        
        // Start updating label positions
        requestAnimationFrame(updateLabelPositions);
    }

    function formatRegionName(name) {
        return name.replace(/_/g, ' ').replace(/\b\w/g, function(l) { return l.toUpperCase(); });
    }

    function updateLabelPositions() {
        if (!window.viewer || !showLabels || !window.viewer.camera) {
            requestAnimationFrame(updateLabelPositions);
            return;
        }
        
        var THREE = BrainBrowser.SurfaceViewer.THREE;
        var container = document.getElementById('brainbrowser');
        if (!container) {
            requestAnimationFrame(updateLabelPositions);
            return;
        }
        
        var rect = container.getBoundingClientRect();
        var camera = window.viewer.camera;
        
        KEY_REGIONS.forEach(function(regionName) {
            var coords = REGION_COORDS[regionName];
            if (!coords) return;
            
            var label = document.getElementById('label-' + regionName);
            if (!label) return;
            
            try {
                // Project 3D position to 2D screen
                var vector = new THREE.Vector3(coords.x, coords.y, coords.z);
                vector.project(camera);
                
                // Check if in front of camera
                if (vector.z > 1) {
                    label.style.display = 'none';
                    return;
                }
                
                // Convert to screen coordinates
                var x = (vector.x * 0.5 + 0.5) * rect.width;
                var y = (-vector.y * 0.5 + 0.5) * rect.height;
                
                // Check if on screen
                if (x < -50 || x > rect.width + 50 || y < -50 || y > rect.height + 50) {
                    label.style.display = 'none';
                    return;
                }
                
                label.style.display = 'block';
                label.style.left = x + 'px';
                label.style.top = y + 'px';
            } catch (e) {
                // Silently ignore projection errors
            }
        });
        
        requestAnimationFrame(updateLabelPositions);
    }

    function setupControls(viewer) {
        // Wireframe toggle
        var wireframeBtn = document.getElementById('btn-wireframe');
        if (wireframeBtn) {
            var wireframeOn = false;
            wireframeBtn.addEventListener('click', function() {
                wireframeOn = !wireframeOn;
                viewer.setWireframe(wireframeOn);
                this.classList.toggle('active', wireframeOn);
            });
        }
        
        // Reset view
        var resetBtn = document.getElementById('btn-reset');
        if (resetBtn) {
            resetBtn.addEventListener('click', function() {
                viewer.setView("lateral");
                viewer.zoom = 1;
                resetBrainColors(viewer);
                resetMarkerColors();
                document.getElementById('effect-select').value = '';
                document.getElementById('receptor-select').value = '';
                currentEffect = null;
                currentReceptor = null;
                document.getElementById('region-info').innerHTML = '<p class="placeholder">Select an effect or receptor to explore</p>';
            });
        }
        
        // Toggle labels button
        var labelsBtn = document.getElementById('btn-labels');
        if (labelsBtn) {
            labelsBtn.addEventListener('click', function() {
                showLabels = !showLabels;
                this.classList.toggle('active', showLabels);
                var overlay = document.getElementById('region-labels-overlay');
                if (overlay) overlay.style.display = showLabels ? 'block' : 'none';
                
                // Toggle marker visibility
                regionMarkers.forEach(function(marker) {
                    marker.visible = showLabels;
                });
                viewer.updated = true;
            });
            labelsBtn.classList.add('active'); // Labels on by default
        }
        
        // Effect select
        var effectSelect = document.getElementById('effect-select');
        if (effectSelect) {
            effectSelect.addEventListener('change', function() {
                var effect = this.value;
                if (effect) {
                    // Clear receptor selection when selecting effect
                    document.getElementById('receptor-select').value = '';
                    currentReceptor = null;
                    highlightEffect(viewer, effect);
                } else {
                    currentEffect = null;
                    resetBrainColors(viewer);
                    clearInfoPanel();
                }
            });
        }
        
        // Receptor select
        var receptorSelect = document.getElementById('receptor-select');
        if (receptorSelect) {
            receptorSelect.addEventListener('change', function() {
                var receptor = this.value;
                if (receptor) {
                    // Clear effect selection when selecting receptor
                    document.getElementById('effect-select').value = '';
                    currentEffect = null;
                    highlightReceptorDensity(viewer, receptor);
                } else {
                    currentReceptor = null;
                    resetBrainColors(viewer);
                    clearInfoPanel();
                }
            });
        }
    }

    function populateReceptorSelect() {
        var select = document.getElementById('receptor-select');
        console.log("populateReceptorSelect called, select:", select, "DATA:", window.DATA);
        if (!select || !window.DATA) {
            console.log("Missing select or DATA, retrying in 500ms");
            setTimeout(populateReceptorSelect, 500);
            return;
        }
        
        var receptors = window.DATA.receptors;
        console.log("Receptors to populate:", Object.keys(receptors).length);
        
        while (select.options.length > 1) {
            select.remove(1);
        }
        
        Object.keys(receptors).sort().forEach(function(receptorName) {
            var option = document.createElement('option');
            option.value = receptorName;
            option.textContent = receptorName;
            select.appendChild(option);
        });
        console.log("Receptor select populated with", select.options.length - 1, "options");
    }

    var effectSelectPopulated = false;
    
    function populateEffectSelect() {
        var select = document.getElementById('effect-select');
        if (!select || !window.DATA) {
            setTimeout(populateEffectSelect, 500);
            return;
        }
        
        // Skip if already populated
        if (effectSelectPopulated && select.options.length > 1) {
            return;
        }
        
        var effects = window.DATA.effects;
        if (!effects) return;
        
        // Clear all existing options except first
        while (select.options.length > 1) {
            select.remove(1);
        }
        
        // Helper to check if an effect has any associated drugs
        function effectHasDrugs(effectKey) {
            if (typeof EFFECT_DRUG_RELATIONSHIPS === 'undefined') return true;
            
            // Direct match
            if (EFFECT_DRUG_RELATIONSHIPS[effectKey] && EFFECT_DRUG_RELATIONSHIPS[effectKey].length > 0) {
                return true;
            }
            
            // Try alternate key formats
            var altKey = effectKey.replace(/_/g, '-');
            if (EFFECT_DRUG_RELATIONSHIPS[altKey] && EFFECT_DRUG_RELATIONSHIPS[altKey].length > 0) {
                return true;
            }
            
            // Try partial match
            var effectLower = effectKey.toLowerCase();
            var keys = Object.keys(EFFECT_DRUG_RELATIONSHIPS);
            for (var i = 0; i < keys.length; i++) {
                if (keys[i].toLowerCase().indexOf(effectLower) !== -1 || 
                    effectLower.indexOf(keys[i].toLowerCase()) !== -1) {
                    if (EFFECT_DRUG_RELATIONSHIPS[keys[i]].length > 0) {
                        return true;
                    }
                }
            }
            
            return false;
        }
        
        // Define PsychonautWiki category order
        var categoryOrder = [
            'Visual effects',
            'Auditory effects', 
            'Tactile effects',
            'Disconnective effects',
            'Smell & taste effects',
            'Multisensory effects',
            'Cognitive effects',
            'Physical effects',
            'Uncomfortable physical effects'
        ];
        
        // Group effects by category (deduplicated by key)
        var categories = {};
        var seenKeys = {};
        var skippedNoDrugs = 0;
        Object.keys(effects).forEach(function(effectKey) {
            var effect = effects[effectKey];
            if (!effect || !effect.name) return;
            
            // Skip duplicates
            if (seenKeys[effectKey]) return;
            seenKeys[effectKey] = true;
            
            // Skip effects without any associated drugs
            if (!effectHasDrugs(effectKey)) {
                skippedNoDrugs++;
                return;
            }
            
            // Normalize category: trim and collapse whitespace
            var category = (effect.category || 'Other').replace(/\s+/g, ' ').trim();
            if (!categories[category]) categories[category] = [];
            
            categories[category].push({
                key: effectKey,
                name: effect.name.trim()
            });
        });
        
        // Sort effects within each category and remove duplicates
        Object.keys(categories).forEach(function(cat) {
            // Deduplicate by key
            var seen = {};
            categories[cat] = categories[cat].filter(function(e) {
                if (seen[e.key]) return false;
                seen[e.key] = true;
                return true;
            });
            // Sort
            categories[cat].sort(function(a, b) {
                return a.name.localeCompare(b.name);
            });
        });
        
        // Add optgroups in order
        categoryOrder.forEach(function(category) {
            if (!categories[category] || categories[category].length === 0) {
                console.log("Skipping empty category:", category);
                return;
            }
            
            console.log("Adding category:", category, "with", categories[category].length, "effects");
            var optgroup = document.createElement('optgroup');
            optgroup.label = category;
            
            categories[category].forEach(function(effect) {
                var option = document.createElement('option');
                option.value = effect.key;
                option.textContent = effect.name;
                optgroup.appendChild(option);
            });
            
            select.appendChild(optgroup);
        });
        
        // Add any uncategorized effects
        Object.keys(categories).forEach(function(category) {
            if (categoryOrder.indexOf(category) === -1 && categories[category].length > 0) {
                var optgroup = document.createElement('optgroup');
                optgroup.label = category;
                
                categories[category].forEach(function(effect) {
                    var option = document.createElement('option');
                    option.value = effect.key;
                    option.textContent = effect.name;
                    optgroup.appendChild(option);
                });
                
                select.appendChild(optgroup);
            }
        });
        
        effectSelectPopulated = true;
        var totalShown = Object.keys(effects).length - skippedNoDrugs;
        console.log("Effect select populated with", totalShown, "effects (filtered", skippedNoDrugs, "with no drugs)");
    }

    function populateDrugSelect() {
        // Could add a drug selector dropdown if needed
    }

    // Pre-compute global average binding profile across all drugs
    var globalAverageBinding = null;
    
    function computeGlobalAverageBinding() {
        if (!window.DATA || !window.DATA.drugs) return {};
        
        var receptorSums = {};
        var receptorCounts = {};
        
        Object.keys(window.DATA.drugs).forEach(function(drugName) {
            var drug = window.DATA.drugs[drugName];
            if (drug.receptors) {
                Object.keys(drug.receptors).forEach(function(receptor) {
                    if (!receptorSums[receptor]) {
                        receptorSums[receptor] = 0;
                        receptorCounts[receptor] = 0;
                    }
                    receptorSums[receptor] += drug.receptors[receptor];
                    receptorCounts[receptor]++;
                });
            }
        });
        
        var averages = {};
        Object.keys(receptorSums).forEach(function(receptor) {
            averages[receptor] = receptorSums[receptor] / receptorCounts[receptor];
        });
        
        return averages;
    }

    function highlightEffect(viewer, effectName) {
        if (!window.DATA) return;
        
        var effect = window.DATA.effects[effectName];
        if (!effect) return;
        
        currentEffect = effectName;
        currentReceptor = null;
        
        // Compute global average if not done yet
        if (!globalAverageBinding) {
            globalAverageBinding = computeGlobalAverageBinding();
        }
        
        // Get drugs that cause this effect
        var drugsForEffect = findDrugsForEffect(effectName);
        
        if (drugsForEffect.length === 0) {
            resetBrainColors(viewer);
            resetMarkerColors();
            updateEffectInfoPanel(effectName, effect, []);
            return;
        }
        
        // Calculate average binding for drugs with this effect
        var effectReceptorSums = {};
        var effectReceptorCounts = {};
        
        drugsForEffect.forEach(function(drugInfo) {
            if (drugInfo.allReceptors) {
                drugInfo.allReceptors.forEach(function(r) {
                    if (!effectReceptorSums[r.name]) {
                        effectReceptorSums[r.name] = 0;
                        effectReceptorCounts[r.name] = 0;
                    }
                    effectReceptorSums[r.name] += r.affinity;
                    effectReceptorCounts[r.name]++;
                });
            }
        });
        
        // Calculate delta from global average (more negative delta = stronger binding than average)
        var receptorDeltas = [];
        Object.keys(effectReceptorSums).forEach(function(receptor) {
            var effectAvg = effectReceptorSums[receptor] / effectReceptorCounts[receptor];
            var globalAvg = globalAverageBinding[receptor] || 0;
            var delta = effectAvg - globalAvg; // More negative = stronger than average
            
            receptorDeltas.push({
                receptor: receptor,
                delta: delta,
                effectAvg: effectAvg,
                globalAvg: globalAvg
            });
        });
        
        // Sort by delta (most negative first = most selective for this effect)
        receptorDeltas.sort(function(a, b) { return a.delta - b.delta; });
        
        // Get brain regions weighted by receptor selectivity (negative delta)
        var regionDensities = {};
        
        receptorDeltas.forEach(function(rd) {
            // Only consider receptors with stronger-than-average binding (negative delta)
            if (rd.delta >= 0) return;
            
            var receptor = window.DATA.receptors[rd.receptor];
            if (receptor && receptor.regions) {
                var weight = Math.abs(rd.delta); // Use magnitude of delta as weight
                
                receptor.regions.forEach(function(regionData) {
                    var regionName = regionData.region;
                    if (!regionDensities[regionName]) {
                        regionDensities[regionName] = 0;
                    }
                    regionDensities[regionName] += regionData.density * weight;
                });
            }
        });
        
        // Convert to array for coloring
        var regionsWithDensity = Object.keys(regionDensities).map(function(regionName) {
            return { region: regionName, density: regionDensities[regionName] };
        });
        
        if (regionsWithDensity.length > 0) {
            // Normalize
            var maxDensity = Math.max.apply(Math, regionsWithDensity.map(function(r) { return r.density; }));
            regionsWithDensity.forEach(function(r) {
                r.density = (r.density / maxDensity) * 100;
            });
            
            colorBrainByDensity(viewer, regionsWithDensity);
            highlightMarkers(regionsWithDensity);
        } else {
            resetBrainColors(viewer);
            resetMarkerColors();
        }
        
        // Update info panel with effect details and top selective receptors
        updateEffectInfoPanel(effectName, effect, receptorDeltas.slice(0, 10));
    }

    function findDrugsForEffect(effectName) {
        if (!window.DATA) return [];
        
        // Get drugs from PsychonautWiki/TripSit scraped data only
        if (typeof EFFECT_DRUG_RELATIONSHIPS === 'undefined') {
            return [];
        }
        
        // Convert effect name to lookup key
        var effectKey = effectName.toLowerCase().replace(/ /g, '_').replace(/-/g, '_');
        
        // Also try the raw effect name from the data
        var effect = window.DATA.effects ? window.DATA.effects[effectName] : null;
        var effectDisplayName = effect ? (effect.name || effectName) : effectName;
        var altKey = effectDisplayName.toLowerCase().replace(/ /g, '_').replace(/-/g, '_');
        
        // Look up in scraped data
        var scrapedDrugs = EFFECT_DRUG_RELATIONSHIPS[effectKey] || 
                           EFFECT_DRUG_RELATIONSHIPS[altKey] ||
                           [];
        
        // Also try partial matches
        if (scrapedDrugs.length === 0) {
            Object.keys(EFFECT_DRUG_RELATIONSHIPS).forEach(function(key) {
                if (key.indexOf(effectKey) !== -1 || effectKey.indexOf(key) !== -1) {
                    scrapedDrugs = scrapedDrugs.concat(EFFECT_DRUG_RELATIONSHIPS[key]);
                }
            });
            // Remove duplicates
            scrapedDrugs = scrapedDrugs.filter(function(item, pos) {
                return scrapedDrugs.indexOf(item) === pos;
            });
        }
        
        if (scrapedDrugs.length === 0) {
            return [];
        }
        
        // Return drugs from scraped data with binding info if available in our database
        var result = [];
        scrapedDrugs.forEach(function(drugName) {
            // Find in our drug database for binding info - use strict matching
            var normalizedName = drugName.toLowerCase().replace(/[^a-z0-9]/g, '');
            var matchedDrug = null;
            var matchedDrugName = null;
            
            if (window.DATA.drugs) {
                // First try exact match
                if (window.DATA.drugs[drugName]) {
                    matchedDrug = window.DATA.drugs[drugName];
                    matchedDrugName = drugName;
                } else {
                    // Try normalized exact match only
                    Object.keys(window.DATA.drugs).forEach(function(dbDrugName) {
                        var dbNormalized = dbDrugName.toLowerCase().replace(/[^a-z0-9]/g, '');
                        if (dbNormalized === normalizedName) {
                            matchedDrug = window.DATA.drugs[dbDrugName];
                            matchedDrugName = dbDrugName;
                        }
                    });
                }
            }
            
            var drugClass = 'other';
            var allReceptors = [];
            
            if (matchedDrug) {
                drugClass = matchedDrug.class || 'other';
                if (matchedDrug.receptors) {
                    Object.keys(matchedDrug.receptors).forEach(function(rName) {
                        allReceptors.push({
                            name: rName,
                            affinity: matchedDrug.receptors[rName]
                        });
                    });
                    allReceptors.sort(function(a, b) { return a.affinity - b.affinity; });
                }
            }
            
            // Also check TripSit data for classification
            if (typeof TRIPSIT_DRUG_DATA !== 'undefined') {
                Object.keys(TRIPSIT_DRUG_DATA).forEach(function(tsName) {
                    var tsNorm = tsName.toLowerCase().replace(/[^a-z0-9]/g, '');
                    if (tsNorm === normalizedName) {
                        drugClass = TRIPSIT_DRUG_DATA[tsName].class || drugClass;
                    }
                });
            }
            
            result.push({
                name: drugName,
                class: drugClass,
                allReceptors: allReceptors,
                matchingReceptors: allReceptors.slice(0, 3),
                source: 'PsychonautWiki'
            });
        });
        
        return result;
    }

    function updateEffectInfoPanel(effectName, effect, receptorDeltas) {
        var regionInfo = document.getElementById('region-info');
        if (!regionInfo) return;
        
        // Get proper display name
        var displayName = effect.name || effectName.replace(/_/g, ' ');
        
        // Use source URL from data if available, otherwise construct it
        var sourceUrl = effect.source || ('https://psychonautwiki.org/wiki/' + encodeURIComponent(displayName.replace(/ /g, '_')));
        var sourceLabel = effect.sourceLabel || 'PsychonautWiki';
        
        var html = '<div class="receptor-info effect-info">';
        html += '<h4>✨ ' + displayName + '</h4>';
        
        // Show category breadcrumb
        if (effect.category) {
            var breadcrumb = effect.category;
            if (effect.subcategory) {
                breadcrumb += ' › ' + effect.subcategory;
            }
            html += '<p class="effect-category">' + breadcrumb + '</p>';
        }
        
        // Show description
        if (effect.description && effect.description !== displayName) {
            html += '<p class="receptor-function">' + effect.description + '</p>';
        }
        
        // Source citation with proper link
        html += '<div class="source-citation">';
        html += '<span class="source-label">📚 Source:</span> ';
        html += '<a href="' + sourceUrl + '" target="_blank" class="source-link">' + sourceLabel + '</a>';
        html += '</div>';
        
        // Show selective receptors (delta from average)
        if (receptorDeltas && receptorDeltas.length > 0) {
            var selectiveReceptors = receptorDeltas.filter(function(rd) { return rd.delta < 0; });
            if (selectiveReceptors.length > 0) {
                html += '<div class="receptor-selectivity">';
                html += '<h5>🎯 Receptor Selectivity (vs average drug)</h5>';
                html += '<p class="source-note">Receptors where drugs causing this effect bind stronger than average</p>';
                html += '<div class="selectivity-bars">';
                
                var maxDelta = Math.abs(selectiveReceptors[0].delta);
                selectiveReceptors.slice(0, 8).forEach(function(rd) {
                    var barWidth = Math.min(100, (Math.abs(rd.delta) / maxDelta) * 100);
                    html += '<div class="selectivity-row">';
                    html += '<span class="receptor-name">' + rd.receptor + '</span>';
                    html += '<div class="selectivity-bar-container">';
                    html += '<div class="selectivity-bar" style="width:' + barWidth + '%"></div>';
                    html += '</div>';
                    html += '<span class="delta-value" title="Effect avg: ' + (rd.effectAvg/1000).toFixed(1) + 'k | Global avg: ' + (rd.globalAvg/1000).toFixed(1) + 'k">';
                    html += (rd.effectAvg / 1000).toFixed(1) + 'k / ' + (rd.globalAvg / 1000).toFixed(1) + 'k';
                    html += '</span>';
                    html += '</div>';
                });
                html += '</div></div>';
            }
        }
        
        // Show drugs that produce this effect (from scraped data)
        var drugsWithEffect = findDrugsForEffect(effectName);
        if (drugsWithEffect.length > 0) {
            var isFromPsychonautWiki = drugsWithEffect[0] && drugsWithEffect[0].source === 'PsychonautWiki';
            var drugSectionTitle = '💊 Drugs That Produce This Effect (' + drugsWithEffect.length + ')';
            var drugSourceNote = '<p class="source-note">Drug-effect relationships from <a href="https://psychonautwiki.org" target="_blank">PsychonautWiki</a> + <a href="https://tripsit.me" target="_blank">TripSit</a></p>';
            
            html += '<div class="binding-drugs drugs-detailed">';
            html += '<h5>' + drugSectionTitle + '</h5>';
            html += drugSourceNote;
            html += '<div class="drug-cards">';
            drugsWithEffect.slice(0, 10).forEach(function(drug) {
                var classColor = window.DATA.classColors ? (window.DATA.classColors[drug.class] || '#6b7280') : '#6b7280';
                html += '<div class="drug-card clickable-drug" data-drug="' + drug.name + '" style="cursor:pointer">';
                html += '<div class="drug-card-header">';
                html += '<span class="drug-tag" style="background:' + classColor + '">' + drug.name + '</span>';
                html += '<span class="drug-class">' + drug.class + '</span>';
                html += '</div>';
                
                // Show receptor bindings (highlight effect-related ones)
                if (drug.allReceptors && drug.allReceptors.length > 0) {
                    html += '<div class="drug-receptors">';
                    html += '<span class="receptors-label">Receptors:</span>';
                    drug.allReceptors.slice(0, 6).forEach(function(r) {
                        var isMatch = r.isEffectReceptor;
                        var className = isMatch ? 'receptor-mini receptor-highlight' : 'receptor-mini';
                        html += '<span class="' + className + '">' + r.name + '</span>';
                    });
                    if (drug.allReceptors.length > 6) {
                        html += '<span class="receptor-more">+' + (drug.allReceptors.length - 6) + '</span>';
                    }
                    html += '</div>';
                }
                
                html += '<div class="drug-click-hint">Click for full binding profile →</div>';
                html += '</div>';
            });
            if (drugsWithEffect.length > 10) {
                html += '<div class="more-tag">+' + (drugsWithEffect.length - 10) + ' more compounds</div>';
            }
            html += '</div></div>';
        }
        
        html += '</div>';
        
        regionInfo.innerHTML = html;
        
        // Add click handlers for drug cards
        attachDrugCardClickHandlers(regionInfo);
    }

    function attachDrugCardClickHandlers(container) {
        var drugCards = container.querySelectorAll('.clickable-drug');
        drugCards.forEach(function(card) {
            card.addEventListener('click', function() {
                var drugName = this.getAttribute('data-drug');
                showDrugDockingModal(drugName);
            });
        });
    }

    function showDrugDockingModal(drugName) {
        // Find the drug in our data
        var drugData = null;
        if (window.DATA && window.DATA.drugs) {
            // Try exact match first
            if (window.DATA.drugs[drugName]) {
                drugData = window.DATA.drugs[drugName];
            } else {
                // Try normalized match
                var normalizedSearch = drugName.toLowerCase().replace(/[^a-z0-9]/g, '');
                Object.keys(window.DATA.drugs).forEach(function(dbName) {
                    var dbNormalized = dbName.toLowerCase().replace(/[^a-z0-9]/g, '');
                    if (dbNormalized === normalizedSearch) {
                        drugData = window.DATA.drugs[dbName];
                        drugName = dbName; // Use the actual name from DB
                    }
                });
            }
        }

        // Create modal content
        var modal = document.getElementById('modal');
        var body = document.getElementById('modal-body');
        
        if (!modal || !body) {
            console.error('Modal elements not found');
            return;
        }

        var html = '<h2>' + drugName + '</h2>';
        
        if (drugData) {
            // Drug class badge
            var drugClass = drugData.class || 'other';
            html += '<span class="card-badge ' + drugClass + '" style="display:inline-block;margin-bottom:1rem">' + drugClass + '</span>';
            
            // Full binding profile
            if (drugData.receptors && Object.keys(drugData.receptors).length > 0) {
                html += '<div class="detail-section">';
                html += '<h4>Full Receptor Binding Profile (' + Object.keys(drugData.receptors).length + ' receptors)</h4>';
                html += '<p class="source-note">Binding energies from molecular docking simulations (OpenMM). More negative = stronger binding.</p>';
                
                // Sort by binding energy (most negative first = strongest)
                var sortedReceptors = Object.entries(drugData.receptors)
                    .sort(function(a, b) { return a[1] - b[1]; });
                
                // Find max absolute value for scaling bars
                var maxAbs = Math.max.apply(null, sortedReceptors.map(function(r) { return Math.abs(r[1]); }));
                
                html += '<div class="binding-profile-full">';
                sortedReceptors.forEach(function(entry, index) {
                    var receptor = entry[0];
                    var energy = entry[1];
                    var pct = Math.abs(energy) / maxAbs * 100;
                    var barColor = energy < -80000 ? '#ef4444' : energy < -60000 ? '#f97316' : energy < -40000 ? '#eab308' : '#22c55e';
                    
                    html += '<div class="binding-bar">';
                    html += '<span class="binding-label">' + receptor + '</span>';
                    html += '<div class="binding-value">';
                    html += '<div class="binding-fill" style="width:' + pct + '%;background:' + barColor + '"></div>';
                    html += '</div>';
                    html += '<span class="binding-energy">' + Math.round(energy) + ' kJ/mol</span>';
                    html += '</div>';
                });
                html += '</div>';
            } else {
                html += '<p>No binding data available for this compound.</p>';
            }
        } else {
            html += '<p>Drug data not found in database.</p>';
        }

        // Effects section (from scraped data)
        var drugEffects = [];
        if (typeof EFFECT_DRUG_RELATIONSHIPS !== 'undefined') {
            var searchName = drugName.toLowerCase();
            var searchNormalized = drugName.replace(/[^a-zA-Z0-9]/g, '').toLowerCase();
            
            Object.keys(EFFECT_DRUG_RELATIONSHIPS).forEach(function(effectId) {
                var drugs = EFFECT_DRUG_RELATIONSHIPS[effectId];
                var found = drugs.some(function(d) {
                    var dLower = d.toLowerCase();
                    var dNormalized = d.replace(/[^a-zA-Z0-9]/g, '').toLowerCase();
                    return dLower === searchName || dNormalized === searchNormalized ||
                           dLower.indexOf(searchName) !== -1 || searchName.indexOf(dLower) !== -1;
                });
                
                if (found) {
                    var effectName = effectId.replace(/_/g, ' ');
                    if (typeof EFFECT_METADATA !== 'undefined' && EFFECT_METADATA[effectId]) {
                        effectName = EFFECT_METADATA[effectId].name;
                    }
                    drugEffects.push({ id: effectId, name: effectName });
                }
            });
        }

        if (drugEffects.length > 0) {
            html += '<div class="detail-section">';
            html += '<h4>Subjective Effects (' + drugEffects.length + ')</h4>';
            html += '<p class="source-note">From PsychonautWiki + TripSit</p>';
            html += '<div class="effect-list">';
            drugEffects.slice(0, 25).forEach(function(e) {
                html += '<span class="tag effect-tag">' + e.name + '</span>';
            });
            if (drugEffects.length > 25) {
                html += '<span class="tag more-tag">+' + (drugEffects.length - 25) + ' more</span>';
            }
            html += '</div></div>';
        }

        // Source links
        var pwikiUrl = 'https://psychonautwiki.org/wiki/' + encodeURIComponent(drugName.replace(/ /g, '_'));
        var tripsitUrl = 'https://tripsit.me/factsheets/' + encodeURIComponent(drugName.replace(/ /g, '-').toLowerCase());
        
        html += '<div class="detail-section">';
        html += '<h4>Sources</h4>';
        html += '<div class="source-citations">';
        html += '<div class="source-citation"><span class="source-label">📚</span>';
        html += '<a href="' + pwikiUrl + '" target="_blank">PsychonautWiki - ' + drugName + '</a></div>';
        html += '<div class="source-citation"><span class="source-label">📚</span>';
        html += '<a href="' + tripsitUrl + '" target="_blank">TripSit - ' + drugName + '</a></div>';
        html += '</div></div>';

        body.innerHTML = html;
        modal.classList.add('active');
    }

    function highlightReceptorDensity(viewer, receptorName) {
        if (!window.DATA) return;
        
        var receptor = window.DATA.receptors[receptorName];
        if (!receptor || !receptor.regions) return;
        
        currentReceptor = receptorName;
        
        // Color the brain based on density
        colorBrainByDensity(viewer, receptor.regions);
        
        // Highlight region markers that have density for this receptor
        highlightMarkers(receptor.regions);
        
        // Update info panel with receptor details, drugs, and effects
        updateInfoPanel(receptorName, receptor);
    }

    function highlightMarkers(regions) {
        var THREE = BrainBrowser.SurfaceViewer.THREE;
        var regionNames = regions.map(function(r) { return r.region; });
        var maxDensity = Math.max.apply(Math, regions.map(function(r) { return r.density; }));
        
        regionMarkers.forEach(function(marker) {
            var name = marker.userData.regionName;
            var sphere = marker.children[0];
            var ring = marker.children[1];
            
            // Check if this region has receptor density
            var regionData = null;
            for (var i = 0; i < regions.length; i++) {
                if (regions[i].region === name || 
                    regions[i].region.indexOf(name) !== -1 || 
                    name.indexOf(regions[i].region) !== -1) {
                    regionData = regions[i];
                    break;
                }
            }
            
            if (regionData) {
                var density = regionData.density / maxDensity;
                // High density = red, low = cyan
                var color = new THREE.Color();
                if (density > 0.5) {
                    color.setHSL(0, 1, 0.5); // Red
                } else {
                    color.setHSL(0.5 - density * 0.5, 1, 0.5); // Cyan to yellow
                }
                sphere.material.color = color;
                sphere.material.opacity = 0.9;
                ring.material.color = color;
                ring.material.opacity = 0.5;
                sphere.scale.set(1 + density, 1 + density, 1 + density);
                
                // Highlight label
                var label = document.getElementById('label-' + name);
                if (label) {
                    label.classList.add('active');
                    label.style.borderColor = '#' + color.getHexString();
                }
            } else {
                // Dim regions without this receptor
                sphere.material.color.setHex(0x444444);
                sphere.material.opacity = 0.3;
                ring.material.color.setHex(0x444444);
                ring.material.opacity = 0.1;
                sphere.scale.set(1, 1, 1);
                
                var label = document.getElementById('label-' + name);
                if (label) {
                    label.classList.remove('active');
                    label.style.borderColor = 'rgba(0, 255, 170, 0.3)';
                }
            }
        });
        
        if (window.viewer) window.viewer.updated = true;
    }

    function resetMarkerColors() {
        var THREE = BrainBrowser.SurfaceViewer.THREE;
        
        regionMarkers.forEach(function(marker) {
            var name = marker.userData.regionName;
            var sphere = marker.children[0];
            var ring = marker.children[1];
            
            sphere.material.color.setHex(0x00ffaa);
            sphere.material.opacity = 0.8;
            ring.material.color.setHex(0x00ffaa);
            ring.material.opacity = 0.3;
            sphere.scale.set(1, 1, 1);
            
            var label = document.getElementById('label-' + name);
            if (label) {
                label.classList.remove('active');
                label.style.borderColor = 'rgba(0, 255, 170, 0.3)';
            }
        });
        
        if (window.viewer) window.viewer.updated = true;
    }

    // Map from knowledge graph region names to Destrieux atlas region indices
    var REGION_TO_ATLAS = {
        // Visual cortex areas (V1-V5)
        'v1': [11, 45],  // G_cuneus, S_calcarine (primary visual cortex)
        'v2': [19, 20],  // G_occipital_middle, G_occipital_sup
        'v3': [19, 20, 11],
        'v4': [19, 20, 21],  // G_oc-temp_lat-fusifor
        'v5': [19, 20],  // MT area
        'visual_cortex': [11, 45, 19, 20],  // All visual regions
        'visual_cortex_layer_iv': [11, 45],
        
        // Frontal regions
        'prefrontal_cortex': [15, 16, 12, 13, 14],  // G_front_sup, G_front_middle, G_front_inf-*
        'prefrontal_cortex_layer_ii': [15, 16],
        'prefrontal_cortex_layer_v': [15, 16],
        'orbitofrontal_cortex': [24, 31, 64, 65],  // G_orbital, G_rectus, S_orbital*
        'anterior_cingulate_cortex': [6, 7],  // G_and_S_cingul-Ant, Mid-Ant
        'anterior_cingulate': [6, 7],
        'anterior_cingulate_ba24': [6, 7],
        'frontal_cortex': [15, 16, 12, 13, 14],
        'cortex': [15, 16, 12, 13, 14, 25, 26, 27, 28, 34, 37, 38],  // General cortex
        
        // Motor regions
        'motor_cortex': [29],  // G_precentral
        'premotor_cortex': [29, 69, 70],
        'supplementary_motor_area': [3],  // G_and_S_paracentral
        
        // Temporal regions
        'temporal_cortex': [34, 37, 38, 44],
        'superior_temporal_gyrus': [34, 35, 36],  // G_temp_sup-*
        'inferior_temporal_cortex': [37],  // G_temporal_inf
        'auditory_cortex': [33],  // G_temp_sup-G_T_transv
        
        // Parietal regions
        'parietal_cortex': [25, 26, 27, 28],
        'somatosensory_cortex': [28],  // G_postcentral
        'posterior_parietal_cortex': [27],  // G_parietal_sup
        
        // Occipital regions
        'occipital_cortex': [11, 19, 20, 43],
        
        // Limbic regions
        'hippocampus': [23],  // G_oc-temp_med-Parahip (parahippocampal)
        'hippocampus_ca1': [23],
        'hippocampus_ca1_stratum_oriens': [23],
        'hippocampus_ca1_stratum_radiatum': [23],
        'hippocampus_ca3': [23],
        'hippocampus_dentate': [23],
        'hippocampus_dentate_gyrus': [23],
        'hippocampus_dentate_molecular': [23],
        'hippocampus_pyramidal': [23],
        'hippocampus_pyramidal_layer': [23],
        'amygdala': [23, 44],  // Near temporal pole
        'amygdala_basolateral': [23],
        'amygdala_central': [23],
        'entorhinal_cortex': [23],
        'entorhinal_cortex_layer_ii': [23],
        'cingulate_cortex': [6, 7, 8, 9, 10],  // G_and_S_cingul-*
        'posterior_cingulate': [9, 10],  // G_cingul-Post-*
        'posterior_cingulate_ba23': [9, 10],
        'insula': [17, 18, 48, 49, 50],  // G_Ins_*, S_circular_insula_*
        'insular_cortex': [17, 18, 48, 49, 50],
        'insula_anterior': [17, 18],
        
        // Subcortical (map to nearby cortical areas)
        'thalamus': [42],  // Medial_wall (closest)
        'thalamus_mediodorsal': [42],
        'thalamus_pulvinar': [42],
        'hypothalamus': [42, 32],
        'hypothalamus_paraventricular': [42],
        'striatum': [42],
        'caudate': [42],
        'caudate_nucleus': [42],
        'caudate_striosomes': [42],
        'putamen': [42],
        'nucleus_accumbens': [42, 32],  // G_subcallosal
        'basal_ganglia': [42],
        'globus_pallidus': [42],
        'globus_pallidus_external': [42],
        'globus_pallidus_internal': [42],
        'ventral_tegmental_area': [42],
        'substantia_nigra': [42],
        'substantia_nigra_compacta': [42],
        'habenula': [42],
        'habenula_medial': [42],
        'bed_nucleus_stria_terminalis': [42],
        
        // Brainstem (medial wall)
        'raphe_nuclei': [42],
        'dorsal_raphe': [42],
        'dorsal_raphe_nucleus': [42],
        'median_raphe_nucleus': [42],
        'locus_coeruleus': [42],
        'periaqueductal_gray': [42],
        'brainstem': [42],
        'brainstem_respiratory': [42],
        
        // Cerebellum - not in cortical parcellation, use occipital
        'cerebellum': [43, 2],  // Pole_occipital, G_and_S_occipital_inf
        'cerebellar_cortex': [43, 2],
        'cerebellum_granule': [43, 2],
        'cerebellum_granule_layer': [43, 2],
        'cerebellum_molecular_layer': [43, 2],
        'cerebellum_purkinje': [43, 2],
        
        // Other
        'claustrum': [17, 18],  // Insula nearby
        'septum': [42],
        'septum_lateral': [42],
        'basal_forebrain': [42, 32],
        'olfactory_bulb': [31, 64],  // G_rectus, S_orbital_med-olfact
        'spinal_cord': [42],
        'default_mode_network': [9, 10, 15, 16, 27],  // PCC, mPFC, angular
        'cerebral_cortex_layer_ii_iii': [15, 16, 25, 26],
        'cerebral_cortex_layer_iv': [15, 16, 25, 26],
        'cerebral_cortex_layer_v': [15, 16, 25, 26]
    };

    function applyParcellationColors(viewer) {
        if (!parcellationData) return;
        
        var model = viewer.model;
        if (!model || !model.children || model.children.length === 0) return;
        
        var shape = model.children[0];
        if (!shape || !shape.geometry) return;
        
        var geometry = shape.geometry;
        var positions = geometry.attributes.position.array;
        var numVertices = positions.length / 3;
        
        // Create or get color attribute
        var colors = geometry.attributes.color;
        if (!colors) {
            var colorArray = new Float32Array(numVertices * 4);
            geometry.addAttribute('color', new BrainBrowser.SurfaceViewer.THREE.BufferAttribute(colorArray, 4));
            colors = geometry.attributes.color;
        }
        var colorArray = colors.array;
        
        var labels = parcellationData.vertex_labels;
        var regionColors = parcellationData.region_colors;
        
        // Color each vertex based on its region
        for (var i = 0; i < numVertices && i < labels.length; i++) {
            var label = labels[i];
            var color = regionColors[label] || [128, 128, 128, 255];
            colorArray[i * 4] = color[0] / 255;
            colorArray[i * 4 + 1] = color[1] / 255;
            colorArray[i * 4 + 2] = color[2] / 255;
            colorArray[i * 4 + 3] = 1.0;
        }
        
        colors.needsUpdate = true;
        viewer.updated = true;
    }

    function colorBrainByDensity(viewer, regions) {
        if (!parcellationData) {
            // Fallback to old method if no parcellation
            colorBrainByDensityLegacy(viewer, regions);
            return;
        }
        
        var model = viewer.model;
        if (!model || !model.children || model.children.length === 0) return;
        
        var shape = model.children[0];
        if (!shape || !shape.geometry) return;
        
        var geometry = shape.geometry;
        var colors = geometry.attributes.color;
        if (!colors) return;
        
        var colorArray = colors.array;
        var numVertices = colorArray.length / 4;
        var labels = parcellationData.vertex_labels;
        
        // Find max density for normalization
        var maxDensity = 0;
        regions.forEach(function(r) {
            if (r.density > maxDensity) maxDensity = r.density;
        });
        
        // Helper to normalize region names for lookup
        function normalizeRegionName(name) {
            return name.toLowerCase().replace(/[- ]/g, '_').replace(/_+/g, '_');
        }
        
        // Build a map of atlas region index -> density
        var regionDensityMap = {};
        regions.forEach(function(regionData) {
            var normalizedName = normalizeRegionName(regionData.region);
            var atlasIndices = REGION_TO_ATLAS[normalizedName];
            if (atlasIndices) {
                var normalizedDensity = regionData.density / maxDensity;
                atlasIndices.forEach(function(idx) {
                    if (!regionDensityMap[idx] || regionDensityMap[idx] < normalizedDensity) {
                        regionDensityMap[idx] = normalizedDensity;
                    }
                });
            } else {
                console.log('Unknown region:', regionData.region, '→', normalizedName);
            }
        });
        
        // Color each vertex based on its region's density
        for (var i = 0; i < numVertices && i < labels.length; i++) {
            var label = labels[i];
            var density = regionDensityMap[label];
            
            if (density !== undefined && density > 0.05) {
                // Color gradient: blue -> cyan -> green -> yellow -> red
                var r, g, b;
                if (density < 0.25) {
                    var t = density * 4;
                    r = 0; g = t; b = 1;
                } else if (density < 0.5) {
                    var t = (density - 0.25) * 4;
                    r = 0; g = 1; b = 1 - t;
                } else if (density < 0.75) {
                    var t = (density - 0.5) * 4;
                    r = t; g = 1; b = 0;
                } else {
                    var t = (density - 0.75) * 4;
                    r = 1; g = 1 - t; b = 0;
                }
                colorArray[i * 4] = r;
                colorArray[i * 4 + 1] = g;
                colorArray[i * 4 + 2] = b;
                colorArray[i * 4 + 3] = 1.0;
            } else {
                // Gray for regions without density
                colorArray[i * 4] = 0.4;
                colorArray[i * 4 + 1] = 0.4;
                colorArray[i * 4 + 2] = 0.45;
                colorArray[i * 4 + 3] = 1.0;
            }
        }
        
        colors.needsUpdate = true;
        viewer.updated = true;
    }

    function colorBrainByDensityLegacy(viewer, regions) {
        var model = viewer.model;
        if (!model || !model.children || model.children.length === 0) return;
        
        var shape = model.children[0];
        if (!shape || !shape.geometry) return;
        
        var geometry = shape.geometry;
        var positions = geometry.attributes.position.array;
        var colors = geometry.attributes.color;
        
        if (!colors) {
            // Create color attribute if it doesn't exist
            var colorArray = new Float32Array(positions.length / 3 * 4);
            for (var i = 0; i < colorArray.length; i += 4) {
                colorArray[i] = 0.7;     // R
                colorArray[i + 1] = 0.65; // G
                colorArray[i + 2] = 0.62; // B
                colorArray[i + 3] = 1.0;  // A
            }
            geometry.addAttribute('color', new BrainBrowser.SurfaceViewer.THREE.BufferAttribute(colorArray, 4));
            colors = geometry.attributes.color;
        }
        
        var colorArray = colors.array;
        var numVertices = positions.length / 3;
        
        // Find max density for normalization
        var maxDensity = 0;
        regions.forEach(function(r) {
            if (r.density > maxDensity) maxDensity = r.density;
        });
        
        // Reset to base brain color (grayish)
        for (var i = 0; i < numVertices; i++) {
            colorArray[i * 4] = 0.5;      // R
            colorArray[i * 4 + 1] = 0.5;  // G
            colorArray[i * 4 + 2] = 0.55; // B
            colorArray[i * 4 + 3] = 1.0;  // A
        }
        
        // Color vertices based on region proximity
        for (var i = 0; i < numVertices; i++) {
            var vx = positions[i * 3];
            var vy = positions[i * 3 + 1];
            var vz = positions[i * 3 + 2];
            
            var totalInfluence = 0;
            var weightedR = 0, weightedG = 0, weightedB = 0;
            
            regions.forEach(function(regionData) {
                var regionName = regionData.region;
                var density = regionData.density;
                var normalizedDensity = density / maxDensity;
                
                // Get region coordinates
                var coords = REGION_COORDS[regionName];
                if (!coords) {
                    // Try to find a partial match
                    var keys = Object.keys(REGION_COORDS);
                    for (var k = 0; k < keys.length; k++) {
                        if (regionName.indexOf(keys[k]) !== -1 || keys[k].indexOf(regionName) !== -1) {
                            coords = REGION_COORDS[keys[k]];
                            break;
                        }
                    }
                }
                if (!coords) return;
                
                // Calculate distance from vertex to region center
                var dx = vx - coords.x;
                var dy = vy - coords.y;
                var dz = vz - coords.z;
                var distance = Math.sqrt(dx * dx + dy * dy + dz * dz);
                
                // Gaussian falloff - regions influence nearby vertices
                var sigma = 25; // Spread of influence
                var influence = Math.exp(-(distance * distance) / (2 * sigma * sigma));
                
                if (influence > 0.01) {
                    // Color gradient from blue (low) to red (high) via yellow
                    var r, g, b;
                    if (normalizedDensity < 0.5) {
                        // Blue to cyan to green
                        var t = normalizedDensity * 2;
                        r = 0;
                        g = t;
                        b = 1 - t * 0.5;
                    } else {
                        // Green to yellow to red
                        var t = (normalizedDensity - 0.5) * 2;
                        r = t;
                        g = 1 - t * 0.5;
                        b = 0;
                    }
                    
                    totalInfluence += influence * normalizedDensity;
                    weightedR += r * influence * normalizedDensity;
                    weightedG += g * influence * normalizedDensity;
                    weightedB += b * influence * normalizedDensity;
                }
            });
            
            if (totalInfluence > 0.05) {
                // Blend with base color
                var blend = Math.min(totalInfluence * 2, 1);
                colorArray[i * 4] = colorArray[i * 4] * (1 - blend) + (weightedR / totalInfluence) * blend;
                colorArray[i * 4 + 1] = colorArray[i * 4 + 1] * (1 - blend) + (weightedG / totalInfluence) * blend;
                colorArray[i * 4 + 2] = colorArray[i * 4 + 2] * (1 - blend) + (weightedB / totalInfluence) * blend;
            }
        }
        
        colors.needsUpdate = true;
        viewer.updated = true;
    }

    function resetBrainColors(viewer) {
        // If we have parcellation data, show the region colors
        if (parcellationData) {
            applyParcellationColors(viewer);
            currentReceptor = null;
            return;
        }
        
        var model = viewer.model;
        if (!model || !model.children || model.children.length === 0) return;
        
        var shape = model.children[0];
        if (!shape || !shape.geometry) return;
        
        var colors = shape.geometry.attributes.color;
        if (!colors) return;
        
        var colorArray = colors.array;
        var numVertices = colorArray.length / 4;
        
        // Reset to brain tissue color
        for (var i = 0; i < numVertices; i++) {
            colorArray[i * 4] = 0.9;
            colorArray[i * 4 + 1] = 0.82;
            colorArray[i * 4 + 2] = 0.78;
            colorArray[i * 4 + 3] = 1.0;
        }
        
        colors.needsUpdate = true;
        viewer.updated = true;
        currentReceptor = null;
    }

    function updateInfoPanel(receptorName, receptor) {
        var regionInfo = document.getElementById('region-info');
        if (!regionInfo) return;
        
        var html = '<div class="receptor-info">';
        html += '<h4>' + receptorName + '</h4>';
        html += '<p class="receptor-function">' + (receptor.function || '') + '</p>';
        
        // Add source citations for receptor data
        html += '<div class="source-citations">';
        html += '<div class="source-citation">';
        html += '<span class="source-label">📚 Receptor Data:</span> ';
        html += '<a href="https://www.guidetopharmacology.org/GRAC/ObjectDisplayForward?objectId=' + getGuideToPharmId(receptorName) + '" target="_blank" class="source-link">IUPHAR/BPS Guide to Pharmacology</a>';
        html += '</div>';
        html += '<div class="source-citation">';
        html += '<span class="source-label">🧬 Brain Regions:</span> ';
        html += '<a href="https://www.proteinatlas.org/" target="_blank" class="source-link">Human Protein Atlas</a>';
        html += '</div>';
        html += '</div>';
        
        // Find drugs that bind to this receptor
        var bindingDrugs = findDrugsForReceptor(receptorName);
        if (bindingDrugs.length > 0) {
            html += '<div class="binding-drugs">';
            html += '<h5>🧪 Binding Compounds (' + bindingDrugs.length + ')</h5>';
            html += '<div class="drug-list">';
            bindingDrugs.slice(0, 10).forEach(function(drug) {
                var classColor = window.DATA.classColors ? (window.DATA.classColors[drug.class] || '#6b7280') : '#6b7280';
                var affinityScore = Math.abs(drug.affinity).toFixed(0);
                html += '<div class="drug-item">';
                html += '<span class="drug-tag" style="background:' + classColor + '">' + drug.name + '</span>';
                html += '<span class="affinity-score" title="Docking score (lower = stronger binding)">' + affinityScore + '</span>';
                html += '</div>';
            });
            if (bindingDrugs.length > 10) {
                html += '<div class="more-tag">+' + (bindingDrugs.length - 10) + ' more compounds</div>';
            }
            html += '</div></div>';
        }
        
        // Show effects associated with this receptor
        if (receptor.effects && receptor.effects.length > 0) {
            html += '<div class="receptor-effects">';
            html += '<h5>✨ Associated Effects</h5>';
            html += '<div class="effect-tags">';
            receptor.effects.forEach(function(effect) {
                var effectName = effect.replace(/_/g, ' ');
                html += '<span class="effect-tag clickable-effect" data-effect="' + effect + '">' + effectName + '</span>';
            });
            html += '</div></div>';
        }
        
        // High density regions
        if (receptor.regions && receptor.regions.length > 0) {
            var sortedRegions = receptor.regions
                .slice()
                .sort(function(a, b) { return b.density - a.density; })
                .slice(0, 8);
            
            html += '<div class="density-regions">';
            html += '<h5>🧠 High Density Regions</h5>';
            html += '<ul class="region-list">';
            
            var maxDensity = sortedRegions[0].density;
            sortedRegions.forEach(function(region) {
                var regionName = region.region.replace(/_/g, ' ');
                var barWidth = (region.density / maxDensity * 100);
                html += '<li>';
                html += '<span class="region-name">' + regionName + '</span>';
                html += '<div class="density-bar-container">';
                html += '<div class="density-bar" style="width:' + barWidth + '%"></div>';
                html += '</div>';
                html += '<span class="density-value">' + region.density + '</span>';
                html += '</li>';
            });
            
            html += '</ul></div>';
        }
        
        html += '</div>';
        
        // Color legend
        html += '<div class="color-legend">';
        html += '<h5>Density Scale</h5>';
        html += '<div class="legend-gradient"></div>';
        html += '<div class="legend-labels"><span>Low</span><span>High</span></div>';
        html += '</div>';
        
        regionInfo.innerHTML = html;
    }

    function findDrugsForReceptor(receptorName) {
        if (!window.DATA || !window.DATA.drugs) return [];
        
        var drugs = [];
        Object.keys(window.DATA.drugs).forEach(function(drugName) {
            var drug = window.DATA.drugs[drugName];
            // drug.receptors is an object with receptor names as keys and binding affinities as values
            if (drug.receptors && (receptorName in drug.receptors)) {
                drugs.push({
                    name: drugName,
                    class: drug.class || 'other',
                    affinity: drug.receptors[receptorName] // binding affinity (docking score)
                });
            }
        });
        
        // Sort by binding affinity (more negative = stronger binding)
        drugs.sort(function(a, b) {
            return a.affinity - b.affinity;
        });
        
        return drugs;
    }

    function clearInfoPanel() {
        var regionInfo = document.getElementById('region-info');
        if (regionInfo) {
            regionInfo.innerHTML = '<p class="placeholder">Select a receptor to view density distribution</p>';
        }
    }

})();
