README

This data file is published by the Movebank Data Repository (www.datarepository.movebank.org). As of the time of publication, a version of this published animal tracking dataset can be viewed on Movebank (www.movebank.org) in the study "White-bearded wildebeest in Kenya" (Movebank Study ID 208413731). Individual attributes in the data files are defined below, in the NERC Vocabulary Server at http://vocab.nerc.ac.uk/collection/MVB and in the Movebank Attribute Dictionary at www.movebank.org/cms/movebank-content/movebank-attribute-dictionary. Metadata describing this data package are maintained at https://datacite.org.

This data package includes the following data files:
White-bearded wildebeest in Kenya-accessory.csv
White-bearded wildebeest in Kenya-gps.csv
White-bearded wildebeest in Kenya-reference-data.csv

Data package citation:
Stabach JA, Hughey L, Reid RS, Worden JS, Leimgruber P, Boone RB. 2020. Data from: Comparison of movement strategies of three populations of white-bearded wildebeest. Movebank Data Repository. https://doi.org/10.5441/001/1.h0t27719

These data are described in the following written publications:
Stabach JA, Hughey L, Reid RS, Worden JS, Leimgruber P, Boone RB. In preparation. Comparison of movement strategies of three populations of white-bearded wildebeest.

Stabach JA. 2015. Movement, resource selection, and the physiological stress response of white-bearded wildebeest [dissertation]. [Fort Collins (CO, USA)]: Colorado State University. https://doi.org/10217/167207

-----------

Terms of Use
This data file is licensed by the Creative Commons Zero (CC0 1.0) license. The intent of this license is to facilitate the re-use of works. The Creative Commons Zero license is a "no rights reserved" license that allows copyright holders to opt out of copyright protections automatically extended by copyright and other laws, thus placing works in the public domain with as little legal restriction as possible. However, works published with this license must still be appropriately cited following professional and ethical standards for academic citation.

We highly recommend that you contact the data creator if possible if you will be re-using or re-analyzing data in this file. Researchers will likely be interested in learning about new uses of their data, might also have important insights about how to properly analyze and interpret their data, and/or might have additional data they would be willing to contribute to your project. Feel free to contact us at support@movebank.org if you need assistance contacting data owners.

For additional information, see

License description: http://creativecommons.org/publicdomain/zero/1.0

General Movebank Terms of Use: www.movebank.org/cms/movebank-content/general-movebank-terms-of-use

Movebank user agreement: https://www.movebank.org/cms/movebank-content/data-policy#user_agreement

-----------

Data Attributes
These definitions come from the Movebank Attribute Dictionary, available at http://vocab.nerc.ac.uk/collection/MVB and www.movebank.org/node/2381.

animal death comments: Comments about the death of the animal.
example: hit by a car
units: none
entity described: individual

animal ID: An individual identifier for the animal, provided by the data owner. If the data owner does not provide an Animal ID, an internal Movebank animal identifier is sometimes shown. 
example: 91876A, Gary
units: none
entity described: individual
same as: individual local identifier

animal life stage: The age class or life stage of the animal at the beginning of the deployment. Can be years or months of age or terms such as 'adult', 'subadult' and 'juvenile'. Best practice is to define units in the values if needed (e.g. '2 years').
example: juvenile, adult
units: not defined
entity described: deployment

animal sex: The sex of the animal. Allowed values are 
m = male
f = female
u = unknown
format: controlled list
entity described: individual

animal taxon: The scientific name of the species on which the tag was deployed, as defined by the Integrated Taxonomic Information System (ITIS, www.itis.gov). If the species name can not be provided, this should be the lowest level taxonomic rank that can be determined and that is used in the ITIS taxonomy. Additional information can be provided using the term 'taxon detail'. 
example: Buteo swainsoni
format: controlled list
entity described: individual
same as: individual taxon canonical name

attachment type: The way a tag is attached to an animal. Values are chosen from a controlled list: 
collar = The tag is attached by a collar around the animal's neck
glue = The tag is attached to the animal using glue
harness = The tag is attached to the animal using a harness
implant = The tag is placed under the skin of the animal
tape = The tag is attached to the animal using tape
other = user specified 
format: controlled list
entity described: deployment

deploy off timestamp: The timestamp when the tag deployment ended. 
example: 2009-10-01 12:00:00.000
format: yyyy-MM-dd HH:mm:ss.SSS
units: UTC or GPS time
entity described: deployment
same as: deploy off date

deploy on latitude: The geographic latitude of the location where the animal was released (intended primarily for instances in which the animal release and tag retrieval locations have higher accuracy than those derived from sensor data).
example: 27.3516
units: decimal degrees, WGS84 reference system
entity described: deployment

deploy on longitude: The geographic longitude of the location where the animal was released (intended primarily for instances in which the animal release and tag retrieval locations have higher accuracy than those derived from sensor data).
example: -97.3321
units: decimal degrees, WGS84 reference system
entity described: deployment

deploy on timestamp: The timestamp when the tag deployment started. 
example: 2008-08-30 18:00:00.000 
format: yyyy-MM-dd HH:mm:ss.SSS 
units: UTC or GPS time
entity described: deployment
same as: deploy on date

deployment end comments: A description of the end of a tag deployment, such as cause of mortality or notes on the removal and/or failure of tag.
example: data transmission stopped after 108 days. Cause unknown
units: none
entity described: deployment

deployment end type: A categorical classification of the tag deployment end, from a controlled list:
captured = The tag remained on the animal but the animal was captured or confined
dead = The deployment ended with the death of the animal that was carrying the tag
equipment failure = The tag stopped working
fall off = The attachment of the tag to the animal failed, and it fell of accidentally
other = other
released = The tag remained on the animal but the animal was released from captivity or confinement
removal = The tag was purposefully removed from the animal
unknown = The deployment ended by an unknown cause
format: controlled list
units: none
entity described: deployment

deployment ID: A unique identifier for the deployment of a tag on animal, provided by the data owner. If the data owner does not provide a Deployment ID, an internal Movebank deployment identifier may sometimes be shown. 
example: Jane-Tag42 
units: none 
entity described: deployment

duty cycle: Remarks associated with the duty cycle of a tag during the deployment, describing the times it is on/off and the frequency at which it transmits or records data. Units and time zones should be defined in the remarks.
example: 15-min fixes from 8:00-18:00 local time (0:00-10:00 UTC)
units: not defined
entity described: deployment

event ID: An identifier for the set of values associated with each event, i.e. sensor measurement. A unique event ID is assigned to every time-location or other time-measurement record in Movebank. If multiple measurements are included within a single row of a data file, they will share an event ID. If users import the same sensor measurement to Movebank multiple times, a separate event ID will be assigned to each. 
example: 6340565
units: none
entity described: event

external temperature: The temperature measured by the tag (different from ambient temperature or internal body temperature of the animal).
example: 32.1
units: degrees Celsius
entity described: event
same as: temperature external

GPS DOP: Dilution of precision provided by the GPS.
example: 1.8
units: unitless
entity described: event

GPS fix type raw: The type of GPS fix as provided by the tag. Common values are 1D = no fix
2D = altitude typically not valid
3D = altitude typically valid. Other values might be present. Where possible, use the more standardized 'GPS fix type'.
example: 3D
units: none
entity described: event

height above ellipsoid: The estimated height above the ellipsoid, typically estimated by the tag. If altitudes are calculated as height above mean sea level, use 'height above mean sea level'.
example: 24.8
units: meters
entity described: event

location lat: The geographic latitude of the location as estimated by the sensor. Positive values are east of the Greenwich Meridian, negative values are west of it. example: -41.0982423 
units: decimal degrees, WGS84 reference system 
entity described: event

location long: The geographic longitude of the location as estimated by the sensor. Positive values are east of the Greenwich Meridian, negative values are west of it. 
example: -121.1761111 
units: decimal degrees, WGS84 reference system 
entity described: event

manipulation type: The way in which the animal was manipulated during the deployment. Values are chosen from a controlled list: 
confined = The animal's movement was restricted to within a defined area
none = The animal received no treatment other than the tag attachment
relocated = The animal was released from a site other than the one at which it was captured
manipulated other = The animal was manipulated in some other way, such as a physiological manipulation. 
format: controlled list 
entity described: deployment

sensor type: The type of sensor with which data were collected. All sensors are associated with a tag id, and tags can contain multiple sensor types. Each event record in Movebank is assigned one sensor type. If values from multiple sensors are reported in a single event, the primary sensor is used. Values are chosen from a controlled list: 
acceleration = The sensor collects acceleration data
accessory measurements = The sensor collects accessory measurements, such as battery voltage
Argos Doppler shift = The sensor uses Argos Doppler shift to determine position
barometer = The sensor records air or water pressure
bird ring = The animal is identified by a band or ring that has a unique identifier
GPS = The sensor uses GPS to determine location
magnetometer = The sensor records the magnetic field
natural mark = The animal is identified by a unique natural marking
orientation = Quaternion components describing the orientation of the tag are derived from accelerometer and gyroscope measurements
radio transmitter = The sensor is a classical radio transmitter
solar geolocator = The sensor collects light levels, which are used to determine position (for processed locations)
solar geolocator raw = The sensor collects light levels, which are used to determine position (for raw light-level measurements)
solar geolocator twilight = The sensor collects light levels, which are used to determine position (for twilights calculated from light-level measurements). 
format: controlled list
entity described: event

study name: The name of the study in Movebank. 
example: Coyotes, Kays and Bogan, Albany NY
units: none
entity described: study, event

study site: A location such as the deployment site or colony, or a location-related group such as the herd or pack name.
example: Pickerel Island North
units: none
entity described: deployment

tag ID: A unique identifier for the tag, provided by the data owner. If the data owner does not provide a tag ID, an internal Movebank tag identifier may sometimes be shown. example: 2342 
units: none 
entity described: tag
same as: tag local identifier

tag manufacturer name: The company or person that produced the tag. 
example: Holohil
units: none
entity described: tag

tag model: The model of the tag. 
example: T61 
units: none 
entity described: tag

tag tech spec: Values for a tag-specific technical specification. Can be used to store measurements not described by existing Movebank attributes. Best practice is to define the values and units in the reference data or study details.
example: 8.31
units: not defined
entity described: event
THIS DATASET: Values are activity X and activity Y, separated by ":". Activity values are provided by Lotek collars and represent results summarized in 5-minute intervals.

tag voltage: The voltage as reported by the tag.
example: 2895
units: mV
entity described: event

timestamp: The date and time corresponding to a sensor measurement or an estimate derived from sensor measurements. 
example: 2008-08-14 18:31:00.000 
format: yyyy-MM-dd HH:mm:ss.SSS 
units: UTC or GPS time 
entity described: event

visible: Determines whether an event is visible on the Movebank map. Values are calculated automatically, with TRUE indicating the event has not been flagged as an outlier by 'algorithm marked outlier', 'import marked outlier' or 'manually marked outlier', or that the user has overridden the results of these outlier attributes using 'manually marked valid' = TRUE. Allowed values are TRUE or FALSE. 
units: none 
entity described: event

-----------

More Information
For more information about this repository, see www.movebank.org/cms/movebank-content/data-repository, the FAQ at www.movebank.org/cms/movebank-content/data-repository-faq, or contact us at support@movebank.org.