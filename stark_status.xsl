<?xml version="1.0"?>
<!DOCTYPE xsl:stylesheet [
<!ENTITY nbsp "&#160;">
]>
<!--

    Copyright (C) 2014  Sergey Lamzin, https://github.com/sergeylamzin/stark

    This file is part of the StarK genome assembler.

    StarK is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    StarK is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.

-->
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
	<xsl:variable name="scale" select="200"/>
	<xsl:template match="/stark_status">
		<html>
			<head>
				<style>
					.nowrap {
						white-space: nowrap;
					}
				</style>
			<title>
				<xsl:choose>
					<xsl:when test="PBS_job/Job_Name"><xsl:value-of select="PBS_job/Job_Name" /> - <xsl:value-of select="PBS_job/@id" /></xsl:when>
					<xsl:otherwise>StarK <xsl:value-of select="version"/></xsl:otherwise>
				</xsl:choose>
			</title>
			</head>
			<body>
				<xsl:apply-templates />
			</body>
		</html>
	</xsl:template>
	
	<xsl:template match="version">
		<h2>StarK <xsl:value-of select="."/></h2>
	</xsl:template>
	
	<xsl:template match="stark_coverage_cache_token">
		<h3>Synchronisation Token</h3>
		<table border="1">
			<tr>
				<td>Mask</td>
				<xsl:for-each select="token">
						<td><xsl:value-of select="@id"/></td>
				</xsl:for-each>
			</tr>
			<tr>
				<td><xsl:value-of select="@mask"/></td>
				<xsl:for-each select="token">
						<td><xsl:value-of select="."/></td>
				</xsl:for-each>
			</tr>
		</table>
		<br />
	</xsl:template>
	
	<xsl:template match="stark_threads">
		<table border="1">
			<thead>
				<tr>
					<td>Thread</td>
					<td>Action</td>
					<td>Reads Inserted<xsl:apply-templates select="reads" /></td>
					<td>Avg Insert Time</td>
					<td>CPU</td>
					<td>HADQ</td>
					<xsl:if test="../stark_coverage_cache_token"><td>Token</td></xsl:if>
				</tr>
			</thead>
			<xsl:apply-templates select="stark_thread_status" />
		</table>
	</xsl:template>
	
	<xsl:template match="reads">
		<br />
		<xsl:value-of select="@done"/>/<xsl:value-of select="@total"/> (<xsl:value-of select='format-number(100 * @done div @total, "0")'/>%)
	</xsl:template>
	
	<xsl:template match="stark_thread_status">
		<tr>
			<td><xsl:value-of select="@pthread_id"/></td>
			<td class="nowrap"><xsl:value-of select="current_action"/></td>
			<td><xsl:value-of select="reads_inserted"/></td>
			<td><xsl:value-of select="avg_insert_time/ns"/>ns</td>
			<td><xsl:value-of select="proc_status/sched_getcpu"/></td>
			<td><xsl:value-of select="ha_dispatch/size"/></td>
			<xsl:if test="stark_coverage_cache_mytoken"><td><xsl:value-of select="stark_coverage_cache_mytoken/Token"/>, <xsl:value-of select="stark_coverage_cache_mytoken/Index"/></td></xsl:if>
			<xsl:for-each select="proc_status/*">
				<td class="nowrap"><xsl:value-of select="name()"/>:<br /><xsl:value-of select="."/></td>
			</xsl:for-each>
		</tr>
	</xsl:template>

	<xsl:template match="threadpool">
		<h3>Thread Pool</h3>
		<table border="1">
			<tr>
				<td>Num Threads</td>
				<td><xsl:value-of select="num_threads"/></td>
			</tr>
			<tr>
				<td>Num Jobs</td>
				<td><xsl:value-of select="numjobs"/></td>
			</tr>
			<tr>
				<td>Running Threads</td>
				<td><xsl:value-of select="running_threads"/></td>
			</tr>
			<tr>
				<td>Threads:</td>
				<td>
					<table>
						<xsl:for-each select="thread">
							<tr><td><xsl:value-of select="."/></td></tr>
						</xsl:for-each>
					</table>
					
				</td>
			</tr>
		</table>
	</xsl:template>
	
	<xsl:template match="PBS_job">
		<h3>PBS Job</h3>
		<table border="1">
			<xsl:for-each select="*">
				<tr>
					<td>
						<xsl:value-of select="name()"/>
					</td>
					<td>
						<xsl:apply-templates select="." />
					</td>
				</tr>
			</xsl:for-each>
		</table>
	</xsl:template>
	
	<xsl:template match="resources_used">
		<table>
			<xsl:for-each select="*">
				<tr>
					<td>
						<xsl:value-of select="name()"/>
					</td>
					<td>
						<xsl:apply-templates />
					</td>
				</tr>
			</xsl:for-each>
		</table>
	</xsl:template>
	
	<xsl:template match="Resource_List">
		<table>
			<xsl:for-each select="*">
				<tr>
					<td>
						<xsl:value-of select="name()"/>
					</td>
					<td>
						<xsl:apply-templates />
					</td>
				</tr>
			</xsl:for-each>
		</table>
	</xsl:template>
	
	<xsl:template match="proc_status">
	</xsl:template>
	
	<xsl:template match="distributor">
	</xsl:template>

</xsl:stylesheet>
