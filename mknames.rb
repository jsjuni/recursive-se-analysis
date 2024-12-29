#!/usr/bin/env ruby

TYPE = %w{System Segment Element Subsystem Assembly Subassembly Part}
typec = Hash.new { |h, k| h[k] = 0 }

require 'csv'

csv = CSV.new(ARGF)
out = CSV.new(STDOUT)

out << ['id', 'key', 'name']

first = true
csv.each do |row|
  if first
    first = false
    next
  end
  id = row[0]
  key = id.split('.').map { |k| k.gsub(/(\d+)/) { |v| '%05d' % v } }.join('.')
  level = id.split('.').size - 2
  type = TYPE[level]
  typec[type] += 1
  name = "#{type} #{typec[type]}"
  out << [id, key, name]
end
