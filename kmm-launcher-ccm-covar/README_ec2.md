## Set up EC2-related configuration 

+	Register an Amazon Elastic Cloud 2 (EC2) account following the instruction [here](http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/get-set-up-for-amazon-ec2.html). 

+	Prepare an [auth.txt](https://github.com/gifford-lab/GERV/blob/master/kmm-launcher-ccm-covar/example/auth.txt) specification file
	
	**No space or tab allowed before and after colon !**
    
    ```
    realm:us-east-1
    region:us-east-1d
    price:3.0
    rsa_key:/cluster/ec2/starcluster.rsa
    access_key:REDACTED
    secret_key:REDACTED
    keyname:starcluster
    mailaddr:your@email.com
    ```
    
    
    #####Useful options:
 
 	+ `realm`: The job submit realms.
    
    + `region`: The job submit regions. You can check the spot prices of a `c3.8xlarge` and pick a cheap region.

    + `price`: The max bid price. Setting this value to Inf will use on-demand allocation (cannot be killed) but will cost a fixed price of ~$.16 / hr. Use this setting only if there's heavy contention, and you cannot wait. $3 is reasonable. Set too low and your jobs will get killed before completed.
    
    + `rsa_key`: The physical location of the key-pair file for your Amazon EC2 account. This file enables the user to remotely communicate with the EC2 instance created. Checkout [here](http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ec2-key-pairs.html#having-ec2-create-your-key-pair) for a instruction to generate your key-pair file.
    
    + `access_key` and `secret_key`: The credentials for your Amazon EC2 account. Checkout [here](http://docs.aws.amazon.com/AWSSimpleQueueService/latest/SQSGettingStartedGuide/AWSCredentials.html) for instruction.
    
    + `keyname`: The name of the key-pair. From your EC2 [console](https://console.aws.amazon.com/ec2/), go to 'Key Pairs' under 'NETWORK & SECURITY' on the left to find the valid key-pair name matched to your rsa_keys indicated above.
    
    + `mailaddr`: sets the email address that gets emailed at the end of a job. The emails will probably get spam-boxed first, so check spam folder.
    
    #####Optional options:
    
    + `itype`: The instance type: valid alternatives are cc2.8xlarge, to use this you must also change the AMI.
    
    + `ami`: The AMI type: you will want to use the HVM image (`ami-864d84ee`) if you use any other instances like cc2.8xlarge or r3.8xlarge.
