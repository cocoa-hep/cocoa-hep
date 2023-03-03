import { Component, OnInit } from '@angular/core';
import {
  EventDataFormat,
  EventDataImportOption,
  EventDisplayService,
} from 'phoenix-ui-components';
import {
  Configuration,
  PresetView,
  PhoenixMenuNode,
  PhoenixLoader,
  StateManager,
} from 'phoenix-event-display';

import { environment } from '../../../environments/environment';
import eventConfig from '../../../../event-config.json';

// import the downloaded configuration from assets
import phoenixMenuConfig from '../../../assets/files/config/atlas-config.json';

@Component({
  selector: 'app-cocoa',
  templateUrl: './cocoa.component.html',
  styleUrls: ['./cocoa.component.scss'],
})
export class ScdComponent implements OnInit {
  phoenixMenuRoot = new PhoenixMenuNode('Phoenix Menu', 'phoenix-menu');
  eventDataImportOptions: EventDataImportOption[] = [
    EventDataFormat.JSON,
    EventDataFormat.JIVEXML,
    EventDataFormat.ZIP,
  ];
  loaded = false;
  loadingProgress = 0;
  // jsRootLoader: JSRootEventLoader;

  constructor(private eventDisplay: EventDisplayService) {}

  ngOnInit() {
    // this.jsRootLoader = new JSRootEventLoader();

    let defaultEvent: { eventFile: string; eventType: string };
    // Get default event from configuration
    if (environment?.singleEvent) {
      defaultEvent = eventConfig;
    } else {
      defaultEvent = {
        //eventFile: 'assets/files/cocoa/events_ttbar_layercolors.json',
        eventFile: 'assets/files/cocoa/events_Wlep_layercolors.json',
        eventType: 'json',
      };
    }

    // Define the configuration
    const configuration: Configuration = {
      eventDataLoader: new PhoenixLoader(),
      presetViews: [
        new PresetView('Left View', [0, 0, -12000], [0, 0, 0], 'left-cube'),
        new PresetView('Center View', [-500, 12000, 0], [0, 0, 0], 'top-cube'),
        new PresetView('Right View', [0, 0, 12000], [0, 0, 0], 'right-cube'),
      ],
      defaultView: [4000, 0, 4000, 0, 0, 0],
      // Set the phoenix menu to be used (defined above)
      phoenixMenuRoot: this.phoenixMenuRoot,
      // Default event data to fallback to if none given in URL
      // Do not set if there should be no event loaded by default
      defaultEventFile: defaultEvent,
    };

    // Initialize the event display
    this.eventDisplay.init(configuration);

    // Load detector geometries

    // Electromagentic calorimeters
    this.eventDisplay.loadOBJGeometry(
      'assets/geometry/COCOA/ECAL1.obj',
      'ECAL1',
      0x6ede8a,
      'ECAL > ECAL1',
      true,
      true,
      false
    );
    this.eventDisplay.loadOBJGeometry(
      'assets/geometry/COCOA/ECAL2.obj',
      'ECAL2',
      0x92e6a7,
      'ECAL > ECAL2',
      true,
      true,
      false
    );
    this.eventDisplay.loadOBJGeometry(
      'assets/geometry/COCOA/ECAL3.obj',
      'ECAL3',
      0xb7efc5,
      'ECAL > ECAL3',
      true,
      true,
      false
    );

    // Hadronic Calorimeters
    this.eventDisplay.loadOBJGeometry(
      'assets/geometry/COCOA/HCAL1.obj',
      'HCAL1',
      0x8189ff,
      'HCAL > HCAL1',
      true,
      true,
      false
    );
    this.eventDisplay.loadOBJGeometry(
      'assets/geometry/COCOA/HCAL2.obj',
      'HCAL2',
      0x96a2ff,
      'HCAL > HCAL2',
      true,
      true,
      false
    );
    this.eventDisplay.loadOBJGeometry(
      'assets/geometry/COCOA/HCAL3.obj',
      'HCAL3',
      0xaeb8ff,
      'HCAL > HCAL3',
      true,
      true,
      false
    );

    this.eventDisplay
      .getLoadingManager()
      .addProgressListener((progress) => (this.loadingProgress = progress));

    // Load the default configuration
    this.eventDisplay.getLoadingManager().addLoadListenerWithCheck(() => {
      console.log('Loading default configuration.');
      this.loaded = true;

      const urlConfig = this.eventDisplay
        .getURLOptionsManager()
        .getURLOptions()
        .get('config');

      if (!urlConfig) {
        const stateManager = new StateManager();
        stateManager.loadStateFromJSON(phoenixMenuConfig);
      }
    });
  }
}
