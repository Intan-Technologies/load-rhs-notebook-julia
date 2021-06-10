using Printf
using Plots
using StaticArrays

mutable struct StimParameters
    stim_step_size
    charge_recovery_current_limit
    charge_recovery_target_voltage
    amp_settle_mode
    charge_recovery_mode
    StimParameters() = new((0 for _ in 1:length(fieldnames(StimParameters)))...)
end

mutable struct ChannelStruct
    native_channel_name::String
    custom_channel_name::String
    native_order::Int16
    custom_order::Int16
    board_stream::Int16
    chip_channel::Int16
    port_name::String
    port_prefix::String
    port_number::Int16
    electrode_impedance_magnitude::Float32
    electrode_impedance_phase::Float32
end

mutable struct SpikeTriggerStruct
    voltage_trigger_mode::Int16
    voltage_threshold::Int16
    digital_trigger_channel::Int16
    digital_edge_polarity::Int16
end

mutable struct Indices
    amplifier::Int
    board_adc::Int
    board_dac::Int
    board_dig_in::Int
    board_dig_out::Int
    Indices() = new((1 for _ in 1:length(fieldnames(Indices)))...)
end

mutable struct Version
    major::Int
    minor::Int
    Version() = new((0 for _ in 1:length(fieldnames(Version)))...)
    #Version(args...) = new(args...)
end

mutable struct Freq
    dsp_enabled::Int16
    actual_dsp_cutoff_frequency::Float32
    actual_lower_bandwidth::Float32
    actual_lower_settle_bandwidth::Float32
    actual_upper_bandwidth::Float32
    desired_dsp_cutoff_frequency::Float32
    desired_lower_bandwidth::Float32
    desired_lower_settle_bandwidth::Float32
    desired_upper_bandwidth::Float32
    notch_filter_frequency::Float32
    desired_impedance_test_frequency::Float32
    actual_impedance_test_frequency::Float32
    amplifier_sample_rate::Float32
    board_adc_sample_rate::Float32
    board_dig_in_sample_rate::Float32
    Freq() = new((0 for _ in 1:length(fieldnames(Freq)))...)
end

mutable struct Header
    version::Version
    frequency_parameters::Freq
    sample_rate::Float32
    eval_board_mode::Int16
    reference_channel::String
    num_amplifier_channels::Int16
    num_board_adc_channels::Int16
    num_board_dac_channels::Int16
    num_board_dig_in_channels::Int16
    num_board_dig_out_channels::Int16
    amp_settle_mode::Int16
    charge_recovery_mode::Int16
    stim_step_size::Float32
    recovery_current_limit::Float32
    recovery_target_voltage::Float32
    dc_amplifier_data_saved::Int16
    notes
    spike_triggers
    amplifier_channels
    dc_amplifier_channels
    stim_channels
    amp_settle_channels
    charge_recovery_channels
    compliance_limit_channels
    board_adc_channels
    board_dac_channels
    board_dig_in_channels
    board_dig_out_channels
    Header() = new(Version(), Freq(), 0, 0, "", 0, 0, 0, 0, 0)
end

mutable struct Data
    t
    amplifier_data
    dc_amplifier_data
    stim_data
    stim_data_raw
    amp_settle_data
    charge_recovery_data
    compliance_limit_data
    stim_polarity
    board_adc_data
    board_dac_data
    board_dig_in_data
    board_dig_in_raw
    board_dig_out_data
    board_dig_out_raw
    Data() = new()
end

mutable struct Result
    t
    spike_triggers
    notes
    frequency_parameters
    stim_parameters
    reference_channel
    amplifier_channels
    amplifier_data
    dc_amplifier_channels
    dc_amplifier_data_saved
    dc_amplifier_data
    stim_channels
    stim_data
    amp_settle_channels
    amp_settle_data
    charge_recovery_channels
    charge_recovery_data
    compliance_limit_channels
    compliance_limit_data
    board_adc_channels
    board_adc_data
    board_dac_channels
    board_dac_data
    board_dig_in_channels
    board_dig_in_data
    board_dig_out_channels
    board_dig_out_data
    Result() = new()
end


function readQString(fid)
    #= Read Qt style String. The first 32-bit unsigned number indicates the length of the string (in bytes).
    If this number equals 0xffffffff, the string is null =#
    
    a = ""
    length = read(fid, UInt32)
    if length == 0xffffffff
        return
    end
    
    # Convert length from bytes to 16-bit Unicode words
    length = length / 2
    for i = 1:length
        thisChar = Char(read(fid, UInt16))
        a = a * thisChar
    end
    return a
end


function plural(n)
    # s = plural(n)
    # Utility function to optionally pluralize words based on the value of n
    if n == 1
        s = ""
    else
        s = "s"
    end
    return s
end


function get_bytes_per_data_block(header)
    # Calculates the number of bytes in each 128-sample datablock.
    N = 128 # n of amplifier samples
    
    # Each data block contains N amplifier samples
    bytes_per_block = N * 4 # timestamp data
    
    bytes_per_block = bytes_per_block + N * 2 * header.num_amplifier_channels
    
    # DC amplifier voltage (absent if flag was off)
    if header.dc_amplifier_data_saved > 0
        bytes_per_block = bytes_per_block + N * 2 * header.num_amplifier_channels
    end
    
    # Stimulation data, one per enabled amplifier channel
    bytes_per_block = bytes_per_block + N * 2 * header.num_amplifier_channels
    
    # Board analog inputs are sampled at same rate as amplifiers
    bytes_per_block = bytes_per_block + N * 2 * header.num_board_adc_channels
    
    # Board analog outputs are sampled at same rate as amplifiers
    bytes_per_block = bytes_per_block + N * 2 * header.num_board_dac_channels
    
    # Board digital inputs are sampled at same rate as amplifiers
    if header.num_board_dig_in_channels > 0
        bytes_per_block = bytes_per_block + N * 2
    end
        
    # Board digital outputs are sampled at same rate as amplifiers
    if header.num_board_dig_out_channels > 0
        bytes_per_block = bytes_per_block + N * 2
    end
    
    return bytes_per_block
end


# Define read_header function
function read_header(fid)
    # Check 'magic number' at beginning of file to make sur ethis is an Intan Technologies RHS2000 data file.
    magic_number = read(fid, UInt32)
    if magic_number != 0xd69127ac
        error("Unrecognized file type.")
    end
    
    header = Header()
    # Read version number
    version = Version()
    version.major = read(fid, Int16)
    version.minor = read(fid, Int16)
    
    println("\nReading Intan Technologies RHS2000 Data File, Version ", version.major, ".", version.minor)
    
    # Read information of sampling rate and amplifier frequency settings
    header.sample_rate = read(fid, Float32)
    header.frequency_parameters.dsp_enabled = read(fid, Int16)
    header.frequency_parameters.actual_dsp_cutoff_frequency = read(fid, Float32)
    header.frequency_parameters.actual_lower_bandwidth = read(fid, Float32)
    header.frequency_parameters.actual_lower_settle_bandwidth = read(fid, Float32)
    header.frequency_parameters.actual_upper_bandwidth = read(fid, Float32)
    header.frequency_parameters.desired_dsp_cutoff_frequency = read(fid, Float32)
    header.frequency_parameters.desired_lower_bandwidth = read(fid, Float32)
    header.frequency_parameters.desired_lower_settle_bandwidth = read(fid, Float32)
    header.frequency_parameters.desired_upper_bandwidth = read(fid, Float32)
    
    # This tells us if a software 50/60 Hz notch filter was enabled during the data acquisition.
    notch_filter_mode = read(fid, Int16)
    notch_filter_frequency = 0
    if notch_filter_mode == 1
        notch_filter_frequency = 50
    elseif notch_filter_mode == 2
        notch_filter_frequency = 60
    end
    header.frequency_parameters.notch_filter_frequency = notch_filter_frequency
    
    header.frequency_parameters.desired_impedance_test_frequency = read(fid, Float32)
    header.frequency_parameters.actual_impedance_test_frequency = read(fid, Float32)
    
    header.amp_settle_mode = read(fid, Int16)
    header.charge_recovery_mode = read(fid, Int16)
    
    header.stim_step_size = read(fid, Float32)
    header.recovery_current_limit = read(fid, Float32)
    header.recovery_target_voltage = read(fid, Float32)
    
    # Place notes in array of Strings
    header.notes = [readQString(fid), readQString(fid), readQString(fid)]
    
    header.dc_amplifier_data_saved = read(fid, Int16)
    header.eval_board_mode = read(fid, Int16)
    
    header.reference_channel = readQString(fid)
    
    # Place frequency-related information in data structure
    header.frequency_parameters.amplifier_sample_rate = header.sample_rate
    
    header.spike_triggers = SpikeTriggerStruct[]
    
    header.amplifier_channels = ChannelStruct[]
    header.dc_amplifier_channels = ChannelStruct[]
    header.stim_channels = ChannelStruct[]
    header.amp_settle_channels = ChannelStruct[]
    header.charge_recovery_channels = ChannelStruct[]
    header.compliance_limit_channels = ChannelStruct[]
    header.board_adc_channels = ChannelStruct[]
    header.board_dac_channels = ChannelStruct[]
    header.board_dig_in_channels = ChannelStruct[]
    header.board_dig_out_channels = ChannelStruct[]
    
    amplifier_index = 1
    board_adc_index = 1
    board_dac_index = 1
    board_dig_in_index = 1
    board_dig_out_index = 1
    
    # Read signal summary from data file header
    number_of_signal_groups= read(fid, Int16)
    
    for signal_group = 1:number_of_signal_groups
        
        signal_group_name = readQString(fid)
        signal_group_prefix = readQString(fid)
        signal_group_enabled = read(fid, Int16)
        signal_group_num_channels = read(fid, Int16)
        signal_group_num_amp_channels = read(fid, Int16)
        
        if (signal_group_num_channels > 0) && (signal_group_enabled > 0)
            for signal_channel = 1:signal_group_num_channels
                new_trigger_channel = SpikeTriggerStruct(0, 0, 0, 0)
                
                new_channel = ChannelStruct("", "", 0, 0, 0, 0, "", "", 0, 0.0, 0.0)
                
                new_channel.port_name = signal_group_name
                new_channel.port_prefix = signal_group_prefix
                new_channel.port_number = signal_group
                
                new_channel.native_channel_name = readQString(fid)
                new_channel.custom_channel_name = readQString(fid)
                new_channel.native_order = read(fid, Int16)
                new_channel.custom_order = read(fid, Int16)
                signal_type = read(fid, Int16)
                channel_enabled = read(fid, Int16)
                new_channel.chip_channel = read(fid, Int16)
                command_stream = read(fid, Int16)
                new_channel.board_stream = read(fid, Int16)
                new_trigger_channel.voltage_trigger_mode = read(fid, Int16)
                new_trigger_channel.voltage_threshold = read(fid, Int16)
                new_trigger_channel.digital_trigger_channel = read(fid, Int16)
                new_trigger_channel.digital_edge_polarity = read(fid, Int16)
                new_channel.electrode_impedance_magnitude = read(fid, Float32)
                new_channel.electrode_impedance_phase = read(fid, Float32)
                
                if channel_enabled > 0
                    if signal_type == 0 
                        push!(header.amplifier_channels, new_channel)
                        
                        # If dc amplifier data is being saved, dc_amplifier_channels
                        if header.dc_amplifier_data_saved > 0
                        new_dc_channel = ChannelStruct("DC_" * new_channel.native_channel_name,
                                "DC_" * new_channel.custom_channel_name,
                                new_channel.native_order,
                                new_channel.custom_order,
                                new_channel.board_stream,
                                new_channel.chip_channel,
                                new_channel.port_name,
                                new_channel.port_prefix,
                                new_channel.port_number,
                                new_channel.electrode_impedance_magnitude,
                                new_channel.electrode_impedance_phase)
                        push!(header.dc_amplifier_channels, new_dc_channel)
                        end
                        
                        # stim_channels
                        new_stim_channel = ChannelStruct("STIM_" * new_channel.native_channel_name,
                            "STIM_" * new_channel.custom_channel_name,
                            new_channel.native_order,
                            new_channel.custom_order,
                            new_channel.board_stream,
                            new_channel.chip_channel,
                            new_channel.port_name,
                            new_channel.port_prefix,
                            new_channel.port_number,
                            new_channel.electrode_impedance_magnitude,
                            new_channel.electrode_impedance_phase)
                        push!(header.stim_channels, new_stim_channel)
                        
                        # amp_settle_channels
                        new_amp_settle_channel = ChannelStruct("AMP_SETTLE_" * new_channel.native_channel_name,
                            "AMP_SETTLE_" * new_channel.custom_channel_name,
                            new_channel.native_order,
                            new_channel.custom_order,
                            new_channel.board_stream,
                            new_channel.chip_channel,
                            new_channel.port_name,
                            new_channel.port_prefix,
                            new_channel.port_number,
                            new_channel.electrode_impedance_magnitude,
                            new_channel.electrode_impedance_phase)
                        push!(header.amp_settle_channels, new_amp_settle_channel)
                        
                        # charge_recovery_channels
                        new_charge_recovery_channel = ChannelStruct("CHARGE_RECOVERY_" * new_channel.native_channel_name,
                            "CHARGE_RECOVERY_" * new_channel.custom_channel_name,
                            new_channel.native_order,
                            new_channel.custom_order,
                            new_channel.board_stream,
                            new_channel.chip_channel,
                            new_channel.port_name,
                            new_channel.port_prefix,
                            new_channel.port_number,
                            new_channel.electrode_impedance_magnitude,
                            new_channel.electrode_impedance_phase)
                        push!(header.charge_recovery_channels, new_charge_recovery_channel)
                        
                        # compliance_limit_channels
                        new_compliance_limit_channel = ChannelStruct("COMPLIANCE_LIMIT_" * new_channel.native_channel_name,
                            "COMPLIANCE_LIMIT_" * new_channel.custom_channel_name,
                            new_channel.native_order,
                            new_channel.custom_order,
                            new_channel.board_stream,
                            new_channel.chip_channel,
                            new_channel.port_name,
                            new_channel.port_prefix,
                            new_channel.port_number,
                            new_channel.electrode_impedance_magnitude,
                            new_channel.electrode_impedance_phase)
                        push!(header.compliance_limit_channels, new_compliance_limit_channel)
                        
                        push!(header.spike_triggers, new_trigger_channel)
                        amplifier_index = amplifier_index + 1
                    elseif signal_type == 1
                        error("Wrong signal type for the rhs format")
                    elseif signal_type == 2
                        error("Wrong signal type for the rhs format")
                    elseif signal_type == 3
                        push!(header.board_adc_channels, new_channel)
                        board_adc_index = board_adc_index + 1
                    elseif signal_type == 4
                        push!(header.board_dac_channels, new_channel)
                        board_dac_index = board_dac_index + 1
                    elseif signal_type == 5
                        push!(header.board_dig_in_channels, new_channel)
                        board_dig_in_index = board_dig_in_index + 1
                    elseif signal_type == 6
                        push!(header.board_dig_out_channels, new_channel)
                        board_dig_out_index = board_dig_out_index + 1
                    else
                        error("Unknown channel type")
                    end
                    
                end
                
            end
            
        end
        
    end
    
    # Summarize contents of data file
    header.num_amplifier_channels = amplifier_index - 1
    header.num_board_adc_channels = board_adc_index - 1
    header.num_board_dac_channels = board_dac_index - 1
    header.num_board_dig_in_channels = board_dig_in_index - 1
    header.num_board_dig_out_channels = board_dig_out_index - 1
    
    return header
    
end

function data_to_result(header, data, data_present)
    # Moves the header and data (if present) into a common object
    result = Result()
    
    stim_parameters = StimParameters()
    stim_parameters.stim_step_size = header.stim_step_size
    stim_parameters.charge_recovery_current_limit = header.recovery_current_limit
    stim_parameters.charge_recovery_target_voltage = header.recovery_target_voltage
    stim_parameters.amp_settle_mode = header.amp_settle_mode
    stim_parameters.charge_recovery_mode = header.charge_recovery_mode
    result.stim_parameters = stim_parameters
    
    result.stim_channels = header.stim_channels
    result.spike_triggers = header.spike_triggers
    result.notes = header.notes
    result.frequency_parameters = header.frequency_parameters
    
    result.reference_channel = header.reference_channel
    result.amplifier_channels = header.amplifier_channels
    result.board_adc_channels = header.board_adc_channels
    result.board_dac_channels = header.board_dac_channels
    
    result.dc_amplifier_data_saved = header.dc_amplifier_data_saved
    if header.dc_amplifier_data_saved > 0
        result.dc_amplifier_channels = header.dc_amplifier_channels
    end
    
    result.compliance_limit_channels = header.compliance_limit_channels
    result.charge_recovery_channels = header.charge_recovery_channels
    result.amp_settle_channels = header.amp_settle_channels
    result.board_dig_in_channels = header.board_dig_in_channels
    result.board_dig_out_channels = header.board_dig_out_channels
    
    if data_present > 0
        result.t = data.t
        result.stim_data = data.stim_data
        result.amplifier_data = data.amplifier_data
        result.board_adc_data = data.board_adc_data
        result.board_dac_data = data.board_dac_data
        if header.dc_amplifier_data_saved > 0
            result.dc_amplifier_data = data.dc_amplifier_data
        end
        result.compliance_limit_data = data.compliance_limit_data
        result.charge_recovery_data = data.charge_recovery_data
        result.amp_settle_data = data.amp_settle_data
        result.board_dig_in_data = data.board_dig_in_data
        result.board_dig_out_data = data.board_dig_out_data
    end
    
    return result
end



function read_one_data_block(data, header, indices, fid)
    # Reads one 128 sample data block from fid into data, at the location indicated by indices
    
    data.t[indices.amplifier:(indices.amplifier + 128 - 1)] = reinterpret(Int32, read(fid, 128 * 4))
    
    if header.num_amplifier_channels > 0
        data.amplifier_data[:, indices.amplifier:(indices.amplifier + 128 - 1)] = permutedims(reshape(reinterpret(UInt16, read(fid, 128 * header.num_amplifier_channels * 2)), Int64(128), Int64(header.num_amplifier_channels)))
        
        # Check if dc amplifier voltage was saved
        if header.dc_amplifier_data_saved > 0
            data.dc_amplifier_data[:, indices.amplifier:(indices.amplifier + 128 - 1)] = permutedims(reshape(reinterpret(UInt16, read(fid, 128 * header.num_amplifier_channels * 2)), Int64(128), Int64(header.num_amplifier_channels)))            
        end
        
        # Get the stimulation data
        data.stim_data_raw[:, indices.amplifier:(indices.amplifier + 128 - 1)] = permutedims(reshape(reinterpret(UInt16, read(fid, 128 * header.num_amplifier_channels * 2)), Int64(128), Int64(header.num_amplifier_channels)))
        
    end
    
    if header.num_board_adc_channels > 0
        data.board_adc_data[:, indices.board_adc:(indices.board_adc + 128 - 1)] = permutedims(reshape(reinterpret(UInt16, read(fid, 128 * header.num_board_adc_channels * 2)), Int64(128), Int64(header.num_board_adc_channels)))
    end
    
    if header.num_board_dac_channels > 0
        data.board_dac_data[:, indices.board_dac:(indices.board_dac + 128 - 1)] = permutedims(reshape(reinterpret(UInt16, read(fid, 128 * header.num_board_dac_channels * 2)), Int64(128), Int64(header.num_board_dac_channels)))
    end
    
    if header.num_board_dig_in_channels > 0
        data.board_dig_in_raw[indices.board_dig_in:(indices.board_dig_in + 128 - 1)] = reinterpret(UInt16, read(fid, 128 * 2))
    end
    
    if header.num_board_dig_out_channels > 0
        data.board_dig_out_raw[indices.board_dig_out:(indices.board_dig_out + 128 - 1)] = reinterpret(UInt16, read(fid, 128 * 2))
    end
end

function notch_filter(input, f_sample, f_notch, bandwidth)
    t_step = 1 / Float64(f_sample)
    f_c = f_notch * t_step

    l = length(input)

    # Calculate IIR filter parameters
    d = exp(-2 * pi * (bandwidth / 2) * t_step)
    b = (1 + d * d) * cos(2 * pi * f_c)
    a0 = 1
    a1 = -b
    a2 = d * d
    a = (1 + d * d) / 2
    b0 = 1
    b1 = -2 * cos(2 * pi * f_c)
    b2 = 1

    output = Vector{Float64}(undef, length(input))
    output[1] = input[1]
    output[2] = input[2]

    #= (If filtering a continuous data stream, change output[1] and output[2] to the previous final two values of out.) =#
    #= Run filter =#
    for k = 3:l
        output[k] = (a*b2*input[k-2] + a*b1*input[k-1] + a*b0*input[k] - a2*output[k-2] - a1*output[k-1])/a0
    end

    return output
end


# Define find_channel_in_group function
function find_channel_in_group(channel_name, signal_group)
    for (count, this_channel) in enumerate(signal_group)
        if this_channel.custom_channel_name == channel_name
            return true, count
        end
    end
    return false, 0
end

# Define find_channel_in_header function
function find_channel_in_header(channel_name, header)
    # Look through all present signal groups
    
    # 1. Look through amplifier_channels
    (channel_found, channel_index) = find_channel_in_group(channel_name, header.amplifier_channels)
    if channel_found
        return true, "amplifier_channels", channel_index
    end
    
    # 2. Look through dc_amplifier_channels
    if header.dc_amplifier_data_saved > 0
        (channel_found, channel_index) = find_channel_in_group(channel_name, header.dc_amplifier_channels)
        if channel_found
            return true, "dc_amplifier_channels", channel_index
        end
    end
    
    # 3. Look through stim_channels
    (channel_found, channel_index) = find_channel_in_group(channel_name, header.stim_channels)
    if channel_found
        return true, "stim_channels", channel_index
    end
    
    # 3.1 Look through amp_settle_channels
    (channel_found, channel_index) = find_channel_in_group(channel_name, header.amp_settle_channels)
    if channel_found
        return true, "amp_settle_channels", channel_index
    end
    
    # 3.2 Look through charge_recovery_channels
    (channel_found, channel_index) = find_channel_in_group(channel_name, header.charge_recovery_channels)
    if channel_found
        return true, "charge_recovery_channels", channel_index
    end
    
    # 3.3 Look through compliance_limit_channels
    (channel_found, channel_index) = find_channel_in_group(channel_name, header.compliance_limit_channels)
    if channel_found
        return true, "compliance_limit_channels", channel_index
    end
    
    # 4. Look through board_adc_channels
    (channel_found, channel_index) = find_channel_in_group(channel_name, header.board_adc_channels)
    if channel_found
        return true, "board_adc_channels", channel_index
    end
    
    # 5. Look through board_dac_channels
    (channel_found, channel_index) = find_channel_in_group(channel_name, header.board_dac_channels)
    if channel_found
        return true, "board_dac_channels", channel_index
    end
    
    # 6. Look through board_dig_in_channels
    (channel_found, channel_index) = find_channel_in_group(channel_name, header.board_dig_in_channels)
    if channel_found
        return true, "board_dig_in_channels", channel_index
    end
    
    # 7. Look through board_dig_out_channels
    (channel_found, channel_index) = find_channel_in_group(channel_name, header.board_dig_out_channels)
    if channel_found
        return true, "board_dig_out_channels", channel_index
    end
    
    return false, "", 0
    
end

# Define plot_channel function
function plot_channel(channel_name, result)
    # Find channel that corresponds to this name
    (channel_found, signal_type, signal_index) = find_channel_in_header(channel_name, result)
    
    # Plot this channel
    if channel_found
        plotly()
        
        if signal_type == "amplifier_channels"
            y_label = "Voltage (microVolts)"
            data_vector = result.amplifier_data[signal_index, :]
            
        elseif signal_type == "dc_amplifier_channels"
            y_label = "Voltage (Volts)"
            data_vector = result.dc_amplifier_data[signal_index, :]
            
        elseif signal_type == "stim_channels"
            y_label = "Current (microAmps)"
            data_vector = result.stim_data[signal_index, :]
            
        elseif signal_type == "amp_settle_channels"
            y_label = "Amp Settle Events (High or Low)"
            data_vector = result.amp_settle_data[signal_index, :]
            
        elseif signal_type == "charge_recovery_channels"
            y_label = "Charge Recovery Events (High or Low)"
            data_vector = result.charge_recovery_data[signal_index, :]
            
        elseif signal_type == "compliance_limit_channels"
            y_label = "Compliance Limit Events (High or Low)"
            data_vector = result.compliance_limit_data[signal_index, :]
            
        elseif signal_type == "board_adc_channels"
            y_label = "Voltage (Volts)"
            data_vector = result.board_adc_data[signal_index, :]
            
        elseif signal_type == "board_dac_channels"
            y_label = "Voltage (Volts)"
            data_vector = result.board_dac_data[signal_index, :]
            
        elseif signal_type == "board_dig_in_channels"
            y_label = "Digital In Events (High or Low)"
            data_vector = result.board_dig_in_data[signal_index, :]
            
        elseif signal_type == "board_dig_out_channels"
            y_label = "Digital Out Events (High or Low)"
            data_vector = result.board_dig_out_data[signal_index, :]
            
        else
            error("Plotting not possible; signal type ", signal_type, " not found")
        end
        
        display(plot(result.t[:], data_vector[:], title = channel_name, xlabel = "Time (s)", ylabel = y_label, legend = false))
        
    else
        error("Plotting not possible; channel ", channel_name, " not found")
    end
end



# Define load_file function
function load_file(filename)
    # Start timing
    start = time()
    
    # Open file
    fid = open(filename, "r")
    filesize = stat(filename).size
    
    # Read file header
    header = read_header(fid)
    
    # Output a summary of recorded data
    println("Found ", header.num_amplifier_channels, " amplifier channel", plural(header.num_amplifier_channels))
    println("Found ", header.num_board_adc_channels, " board ADC channel", plural(header.num_board_adc_channels))
    println("Found ", header.num_board_dac_channels, " board DAC channel", plural(header.num_board_dac_channels))
    println("Found ", header.num_board_dig_in_channels, " board digital input channel", plural(header.num_board_dig_in_channels))
    println("Found ", header.num_board_dig_out_channels, " board digital output channel", plural(header.num_board_dig_out_channels))
    
    # Determine how many samples the data file contains
    bytes_per_block = get_bytes_per_data_block(header)
    
    # Calculate how many data blocks are present
    data_present = 0
    bytes_remaining = filesize - position(fid)
    if bytes_remaining > 0
        data_present = 1
    end
    
    if bytes_remaining % bytes_per_block != 0
        error("Something is wrong with file size: should have a whole number of data blocks")
    end
    
    num_data_blocks = Int(bytes_remaining / bytes_per_block)
    
    # Calculate how many samples of each signal type are present
    num_amplifier_samples = Int(128 * num_data_blocks)
    num_board_adc_samples = Int(128 * num_data_blocks)
    num_board_dac_samples = Int(128 * num_data_blocks)
    num_board_dig_in_samples = Int(128 * num_data_blocks)
    num_board_dig_out_samples = Int(128 * num_data_blocks)
    
    # Calculate how much time has been recorded
    record_time = num_amplifier_samples / header.sample_rate
    
    # Output a summary of contents of header file
    if data_present > 0
        @printf("File contains %0.3f seconds of data. Amplifiers were sampled at %0.2f kS/s.\n", record_time, header.sample_rate / 1000)
    else
        @printf("Header file contains no data. Amplifiers were sampled at %0.2f kS/s.\n", header.sample_rate / 1000)
    end
    
    if data_present > 0
        # Pre-allocate memory for data
        println("Allocating memory for data...\n")
        
        data = Data()
        data.t = zeros(Int32, 1, num_amplifier_samples)
        
        data.amplifier_data = zeros(UInt16, header.num_amplifier_channels, num_amplifier_samples)
        
        if header.dc_amplifier_data_saved > 0
            data.dc_amplifier_data = zeros(UInt16, header.num_amplifier_channels, num_amplifier_samples)
        end
        
        data.stim_data_raw = zeros(UInt16, header.num_amplifier_channels, num_amplifier_samples)
        data.stim_data = zeros(Int16, header.num_amplifier_channels, num_amplifier_samples)
        
        data.board_adc_data = zeros(UInt16, header.num_board_adc_channels, num_board_adc_samples) # Maybe Float64?
        data.board_dac_data = zeros(UInt16, header.num_board_dac_channels, num_board_dac_samples)
        
        # by default, this script interprets digital events (digital inputs and outputs) as Int16
        data.board_dig_in_data = zeros(Int16, header.num_board_dig_in_channels, num_board_dig_in_samples)
        data.board_dig_in_raw = Vector{UInt16}(undef, num_board_dig_in_samples)
        data.board_dig_out_data = zeros(Int16, header.num_board_dig_out_channels, num_board_dig_out_samples)
        data.board_dig_out_raw = Vector{UInt16}(undef, num_board_dig_out_samples)
        
        # Read sampled data from file
        println("Reading data from file...")
        
        # Initialize indices used in looping
        indices = Indices()
        
        print_increment = 10
        percent_done = print_increment
        for i = 1:num_data_blocks
            read_one_data_block(data, header, indices, fid)
            
            # Increment indices
            indices.amplifier = 128 + indices.amplifier
            indices.board_adc = 128 + indices.board_adc
            indices.board_dac = 128 + indices.board_dac
            indices.board_dig_in = 128 + indices.board_dig_in
            indices.board_dig_out = 128 + indices.board_dig_out
            
            fraction_done = 100 * (1.0 * i / num_data_blocks)
            if fraction_done >= percent_done
                println(percent_done, "% done...")
                percent_done = percent_done + print_increment
            end
            
        end
        
        # Make sure we have read exactly the right amount of data
        bytes_remaining = filesize - position(fid)
        if bytes_remaining != 0
            error("Error: End of file not reached.")
        end
    else
        data = nothing
    end
    
    # Close data file
    close(fid)
    
    if data_present > 0
        println("Parsing data...\n")
        
        # Extract digital input channels to separate variables.
        for i = 1 : header.num_board_dig_in_channels
            mask = 2^header.board_dig_in_channels[i].native_order
            data.board_dig_in_data[i,:] = (x -> (x > 0 ? 1 : 0)).(data.board_dig_in_raw .& mask)
        end
        
        # Extract digital output channels to separate variables.
        for i = 1 : header.num_board_dig_out_channels
            mask = 2^header.board_dig_out_channels[i].native_order
            data.board_dig_out_data[i,:] = (x -> (x > 0 ? 1 : 0)).(data.board_dig_out_raw .& mask)
        end
        
        # Extract stimulation data
        data.compliance_limit_data = (x -> (x > 0 ? 1 : 0)).(data.stim_data_raw .& 32768) # get 2^15 bit
        #data.charge_recovery_data = data.stim_data_raw .& 16384 # get 2^14 bit
        data.charge_recovery_data = (x -> (x > 0 ? 1 : 0)).(data.stim_data_raw .& 16384) # get 2^14 bit
        #data.amp_settle_data = data.stim_data_raw .& 8192 # get 2^13 bit
        data.amp_settle_data = (x -> (x > 0 ? 1 : 0)).(data.stim_data_raw .& 8192) # get 2^13 bit
        data.stim_polarity = 1 .- (2 * ((data.stim_data_raw .& 256) .>> 8)) # get 2^8 bit, interpret as +1 for 0_bit or -1 for 1_bit
        
        curr_amp = data.stim_data_raw .& 255 # get least-significant 8 bits corresponding to the current amplitude
        data.stim_data = curr_amp .* data.stim_polarity # multiply current amplitude by the correct sign
        
        # Scale voltage levels appropriately
        data.amplifier_data = 0.195 .* (data.amplifier_data .- 32768) # units = microvolts
        data.stim_data = header.stim_step_size .* (data.stim_data ./ 1.0e-6)
        if header.dc_amplifier_data_saved > 0
            data.dc_amplifier_data = -0.01923 .* (data.dc_amplifier_data .- 512) # units = volts
        end
        data.board_adc_data = 312.5e-6 .* (data.board_adc_data .- 32768) # units = volts
        data.board_dac_data = 312.5e-6 .* (data.board_dac_data .- 32768) # units = volts
        
        # Check for gaps in timestamps.
        num_gaps = sum(diff(data.t, dims=2)[1, :] .!= 1)
        if num_gaps == 0
            println("No missing timestamps in data.")
        else
            println("Warning: ", num_gaps, " gaps in timestamp data found. Time scale will not be uniform!")
        end
        
        # Scale time steps (units = seconds).
        data.t = data.t / header.sample_rate
        
        # If the software notch filter was selected during the recording, apply the same notch filter to amplifier data here
        if header.frequency_parameters.notch_filter_frequency > 0
            println("Applying notch filter...")
            
            print_increment = 10
            percent_done = print_increment
            for i = 1 : header.num_amplifier_channels
                data.amplifier_data[i, :] = notch_filter(data.amplifier_data[i, :], header.sample_rate, header.frequency_parameters.notch_filter_frequency, 10)
                fraction_done = 100 * (i / header.num_amplifier_channels)
                if fraction_done >= percent_done
                    println(percent_done, "% done...")
                    percent_done = percent_done + print_increment
                end
            end
        end
    end
    
    # Move variables to result struct
    result = data_to_result(header, data, data_present)
    
    elapsed = time() - start
    @printf("Done! Elapsed time: %0.1f seconds\n", elapsed)
    
    return(result, data_present)
end

# Define print_all_channel_names function
function print_all_channel_names(result)
    # Print all amplifier_channels
    print_names_in_group(result.amplifier_channels)
    
    # Print all dc_amplifier_channels
    if result.dc_amplifier_data_saved > 0
        print_names_in_group(result.dc_amplifier_channels)
    end
    
    # Print all stim_channels
    print_names_in_group(result.stim_channels)
    
    # Print all amp_settle_channels
    print_names_in_group(result.amp_settle_channels)
    
    # Print all charge_recovery_channels
    print_names_in_group(result.charge_recovery_channels)
    
    # Print all compliance_limit_channels
    print_names_in_group(result.compliance_limit_channels)
    
    # Print all board_adc_channels
    print_names_in_group(result.board_adc_channels)
    
    # Print all board_dac_channels
    print_names_in_group(result.board_dac_channels)
    
    # Print all board_dig_in_channels
    print_names_in_group(result.board_dig_in_channels)
    
    # Print all board_dig_out_channels
    print_names_in_group(result.board_dig_out_channels)
end


# Define function print_names_in_group
function print_names_in_group(signal_group)
    for this_channel in signal_group
        println(this_channel.custom_channel_name)
    end
end