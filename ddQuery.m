%{
                                _     _  ____
                             __| | __| |/ __ \ _   _  ___ _ __ _   _
                            / _` |/ _` | |  | | | | |/ _ \ '__| | | |
                           | (_| | (_| | |__| | |_| |  __/ |  | |_| |
                            \__,_|\__,_|\___\_\\__,_|\___|_|   \__, |
                            easy-as-pie API to Data Dictionary |___/

v0.8, 2017 robert@raschhour.com

ddQuery is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

ddQuery is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along with ddQuery. If not,
see <http://www.gnu.org/licenses/>.
%}

classdef ddQuery < double
	%DDQUERY  easy-as-pie API to Data Dictionary
	
	methods
		function this = ddQuery(query, varargin)
			%if nargin == 0, return, end; % syntax for allocation
			
			if isnumeric(query) % simple handle-array conversion
				handles = query;
			
			elseif iscellstr(query) % pathes-array conversion
				handles = cellfun(@(s) dsdd('GetAttribute', s, 'hDdObject'), query, 'UniformOutput', false);
				%ignore nonmatching: assert(all(~cellfun(@isempty, handles)), 'does not match an object');
				handles = [handles{:}];
				
			elseif ischar(query) % normal query
				% additional arguments represent sets of blocks to pick from via ^123 - index
				idx = cellfun(@isnumeric, varargin);
				varargin(idx) = cellfun(@(hs) {double(reshape(hs, 1, []))}, varargin(idx));
				
				idx = cellfun(@ischar, varargin);
				varargin(idx) = cellfun(@(s) {dsdd('GetAttribute', s, 'hDdObject')}, varargin(idx));
				
				idx = cellfun(@iscellstr, varargin);
				varargin(idx) = cellfun(@(cs) {reshape(cellfun(@(s) dsdd('GetAttribute', s, 'hDdObject'), cs), 1, [])}, varargin(idx));
				
				handles = ddQuery.select(strtrim(query), varargin{:});
			else
				error('illegal arguments');
			end
			
			this = this@double(handles);
		end
	end
	methods (Access='public', Hidden=true)
		function disp(this)
			if isempty(this)
				fprintf('   Empty ddQuery: %d-by-0\n', size(this, 1));
			else
				disp('ddQuery with handles');
				x = dbstack; % hyperlinks don't look nice in datatip displays
				if numel(x) >= 2 && strncmp(x(2).name, 'datatipinfo', 11)
					disp(double(this));
				else
					% must build one big array of text, so that it can print instantaneously
					repr = arrayfun(@(h) sprintf('%*s<a href="matlab: dsddman(''Select'', %d);">%d</a>', 9-floor(log10(h)), '', h, h), ...
						double(this)', 'UniformOutput', false);
%                     repr(end+1, :) = {char(10)}; disp([repr{:}])
					
					w = get(0, 'CommandWindowSize'); w = floor((w(1)-3) / 10);
					if size(repr, 1) > w
						% must print groups of columns that fit the screen
						for i = 1:w:size(repr, 1)
							j = min(i+w-1, size(repr, 1));
							fprintf('  Columns %d through %d\n', i, j)
							crepr = [repr(i:j, :)];
							crepr(end+1, 1:end-1) = {char(10)}; %#ok<AGROW>
							disp([crepr{:}]);
						end
					else % can print evrything in one go.
						repr(end+1, 1:end-1) = {char(10)};
						disp([repr{:}])
					end
				end
			end
		end
		function sel = subsref(this, subs)
			sel = this; % initial selection of dd-nodes
			for sub = subs
				switch sub.type
					case '()' % selection
						if isa(sel, 'ddQuery')
							if isscalar(sub.subs) && isnumeric(sub.subs{1}) && ~isvector(sel)
								sub.subs = [sub.subs, ':'];
							elseif isscalar(sub.subs) && islogical(sub.subs{1})
								sub.subs = [':', sub.subs];
							end
							sel = ddQuery(builtin('subsref', double(sel), sub));
						else
							sel = builtin('subsref', sel, sub);
						end
						
					case '.' % parameter access
						switch lower(sub.subs)
							case 'hyperlink' % an array of links to the dd-objects (with their name)
								sel = arrayfun( ...
									@(h) {sprintf('<a href="matlab: dsddman(''Select'', %d);">%s</a>' ...
									, h, dsdd('GetAttribute', h, 'Name'))}, ...
									double(sel));
								
							case 'handle' % most functions prefer plain double over objects derived from double
								sel = double(sel);
								
							case 'wrap' % TODO: is this cool?
								sel = ddQuery(sel);
								
							case 'parent'
								sel = arrayfun(@(h) {dsdd('GetAttribute', h, 'hDdParent')}, double(sel));
								
							case 'path' % return the path of all blocks
								sel = arrayfun(@(h) {dsdd('GetAttribute', h, 'NormalizedPath')}, double(sel));
								
							case 'name'
								sel = arrayfun(@(h) dsdd('GetAttribute', h, 'Name'), double(sel), 'UniformOutput', false);
								
							case {'kind', 'objectkind'}
								sel = arrayfun(@(h) dsdd('GetAttribute', h, 'ObjectKind'), double(sel), 'UniformOutput', false);
								
							otherwise % select an object property, attribute or field of struct
								if isnumeric(sel) % property of a dd-object
									%arpids = dsdd('GetAutoRenamePropertyIndices', double(sel), sel.subs);
									if isempty(regexp(sub.subs, 'Target$', 'once'))
										sel = arrayfun(@(h) dsdd('Get', h, sub.subs), double(sel), 'UniformOutput', false);
									else
										sel = arrayfun(@(h) dsdd(['Get' sub.subs], h), double(sel), 'UniformOutput', false);
									end
									
								elseif isstruct(sel) % field of a struct
									sel = arrayfun(@(h) h.(sub.subs), sel, 'UniformOutput', false);
									
								end
						end
						
						% TODO: when the result is a cell but happens to be a uniform (e.g.
						% LineHandles.Outport of a bunch of inport blocks), the result cab be a
						% uniform double-array instead of a cell nesting all the single entries.
						% This is counter-intuitive for LineHandles (of arbitrary blocks) but
						% not for e.g. Position
						
					case '{}' % actions
						
					otherwise % they did something stupid
						error(['what is ''' sub.type '''?']);
				end
				
				% coerce the result to a common convenient data type
				if iscell(sel)
					if isscalar(sel) % TODO: this isn't actually cool, because it breaks algorhitms that work with entire query-results, when those just happen to be scalar
						sel = sel{1};
					else
						type = unique(cellfun(@class, sel, 'UniformOutput', false));
						if ~isscalar(type), continue, end

						if ~all(cellfun(@isscalar, sel)), continue, end
						
						if isnumeric(sel{1})
							sel = cell2mat(sel);
						elseif isstruct(sel{1})
							fields = cellfun(@fieldnames, sel, 'UniformOutput', false);
							if ~isequal(fields{:}), return, end
							sel = cell2mat(sel);
						end
					end
				end
			end
		end
		function this = subsasgn(this, subs, value)
			% initial selection of dd-nodes
			sel = this;
			for sub = subs
				switch sub.type
					case '()' % slicing
						if isa(sel, 'ddQuery')
							if isscalar(sub.subs) && isnumeric(sub.subs{1}) && ~isvector(sel)
								sub.subs = [sub.subs, ':'];
							elseif isscalar(sub.subs) && islogical(sub.subs{1})
								sub.subs = [':', sub.subs];
							end
							sel = slQuery(builtin('subsref', double(sel), sub));
						else
							sel = builtin('subsref', sel, sub);
						end
						
					case '.'
						
						if ismember(lower(sub.subs), {'name', 'parent'}) % an attribute, not a property
							
							ddQuery.arrayfun(@(h, v) dsdd('SetAttribute', h, 'name', v), double(sel), value);
						
						else
							
							ddQuery.arrayfun(@(h, v) dsdd('Set', h, sub.subs, v), double(sel), value);
							
						end
					case '{}' % actions
						
					otherwise % they did something stupid
						error(['what is ''' sub.type '''?']);
				end
			end
		end
	end
	
	methods(Access=public) % operators that are operations
		function ps = properties(this)
			ps = sort(dsdd('GetPropertyNames', double(this(1,1))));
		end
		function [varargout] = ctranspose(this) % allow 'dispersed' assignments: [a, b, ~] = ddQuery(...)';
			if nargout < 2
				varargout{1} = this; % this does not transpose the array either
			else
				varargout = cell(1, nargout);
				handles = double(this);
				for i = 1:nargout
					varargout{i} = ddQuery(handles(i, :));
				end
			end
		end
		function new = rdivide(object, spec)
			if ischar(spec) % it's an object-kind specification
				% TODO: parse the spec entirely (param-specs, class)
				spec = regexp(spec, '(?<kind>\w+)(#)?(?<id>(?(2)\w+))', 'names', 'once');
				assert(~isempty(spec));
				if isempty(spec.id), spec.id = ['ddq_' spec.kind]; end
				
				%new = ddQuery(arrayfun(@(o) dsdd(['Add' spec], o, ['ddq_' spec.kind]), double(object)));
				new = ddQuery(ddQuery.arrayfun(@(p) dsdd('Create', [spec.id], 'Parent', p, 'ObjectKind', spec.kind), double(object)));

			elseif isa(spec, 'ddQuery')
				if isscalar(spec)
					new = ddQuery(ddQuery.arrayfun(@(s, p) dsdd('CopyTree', 'Source', s, 'Destination', p), double(spec), double(object)));
				else
					error('not supported')
				end
			end
				
		end
	end
	
	methods(Access=private, Static)
		function handles = select(query, varargin) % core "select" algorithm of ddQuery
			% split query along the combinators                                                                       ( outside [] )
			[selectors, combinators] = regexp(query, '\s*( |\\\\|\\|//|/|(:\s*\w+\s*)?(->|<-|>>|<<)(\s*\w+\s*:)?|,)\s*(?![^\[]*\])', 'split', 'match');
			
			% we start with the combinator '//' for arbitrary descendence and the search root
			
			root = regexprep(dsdd('GetAttribute', dsddman('GetSelected'), 'path'), '^(//DD\d+).*$', '$1'); % root = '//DD0/';
			root = dsdd('GetAttribute', root, 'hddobject'); % always search only in current model
			handles = double.empty(0, 1);
			hot = root;
			for act = [',' combinators; selectors]
				% parse the combinator:     '    (colon with portspec)...(      combinator type (again)                   )...(portspec with colon )
				combinator = regexp(act{1}, '^\s*(:)?(?<sp>(?(1)\w+))?\s*(?<type>( |\\\\|\\|//|/|->|<-|>>|<<|,(?![^\[]*\])))\s*(?<dp>\w+)?\s*(?(4):)?\s*$', 'names');
				
				% cast numeric port qualifiers
				if ~isnan(str2double(combinator.sp)), combinator.sp = str2double(combinator.sp); end
				if ~isnan(str2double(combinator.dp)), combinator.dp = str2double(combinator.dp); end
				
				% build the filter for find from blockspec (all searches are reduced to this)
				find_args = {};
				
				% combinators always represent a search relating to some set of properties from the row of
				% blocks selected previously. In order to avoid multiple seaches based on the same info (that
				% will result in the same result), we can group all rows based on this distinguishing info and
				% perform the search with this and then outer-join the result to the original set of rows
				switch combinator.type
					case ',' % group all into one, because not really a combinator
						hinfos = repmat(root, size(hot));
					case {'\', '\\'} % group by parent of the last object
						% case {'\', '\\', ' '}  NOTE/TODO: the sibling-combinator ' ' can't be here included here because each block cannot be included amongst its own siblings
						hinfos = arrayfun(@(h) dsdd('GetAttribute', h, 'hDdParent'), hot);
					otherwise % group only by object-handle itself
						hinfos = hot;
				end
				
				if regexp(act{2}, '^\$\d+$', 'match', 'once') % selector is a reference ~> simple column-number
					selector = str2double(act{2}(2:end));
					new_handles = double.empty(size(handles, 1), 0); % height stays the same
					
				else % selector is real ~> create structure
					% parse as selector:       ^(*)(    parens around arg index     )?(objectkind)...(   hash with id or special   )...(    brackets and qualifier list   )...( plus and pseudo-class )$
					selector = regexp(act{2}, '^\*?(\()?(?<argidx>(?(1)\d+))?(?(1)\))?(?<kind>\w+)?\s*(#)?(?<id>(?(5)(?:\w+|\.\.?)))?\s*(\[)?(?<attributes>(?(7).+))(?(7)\])\s*(\+)?(?<pseudo>(?(10)\w+))?$', 'names');
					assert(~isempty(selector), 'malformed selector ''%s''', act{2});
					% split the attribute qualifiers:                 '(attribute )...(          operator           )...(    value    )( comma? )
					selector.attributes = regexp(selector.attributes, '(?<name>\w+)\s*(?<operator>(=|\^=|\$=|\*=|~=))\s*(?<value>[^,]+)(\s*,\s*)?', 'names');
					
					if ~isempty(selector.id)
						find_args = [find_args, 'Name', selector.id]; %#ok<AGROW>
					end
					
					if ~isempty(selector.kind)
						find_args = [find_args 'ObjectKind', selector.kind]; %#ok<AGROW>
					end
					
					for attr = selector.attributes % parameters
						if attr.value(1) == '"' && attr.value(end) == '"'
							attr.value = attr.value(2:end-1);
						end
						
						switch attr.operator
							case '=' % literal match "text"
								attr.value = ['^' regexptranslate('escape', attr.value) '$'];
							case '^=' % literal match at beginning "text*"
								attr.value = ['^' regexptranslate('escape', attr.value)];
							case '$=' % literal match at the end "*text"
								attr.value = [regexptranslate('escape', attr.value) '$'];
							case '*=' % literal match anywhere "*text*"
								attr.value = regexptranslate('escape', attr.value);
							case '~=' % regex match
								% attr.value shall be a regexp already
						end
						
						switch lower(attr.name)
							case 'name'
								find_args = [find_args, 'Name', attr.value]; %#ok<AGROW>
							case 'objectkind'
								find_args = [find_args, '$ObjectKindRegexp', attr.value]; %#ok<AGROW> NOTE: this is not a property but magic name for ddQuery.find
							case 'handle'
								error('cannot handle handle');
							otherwise
								%find_args = [find_args, 'Property', {{'name', attr.name, 'value', attr.value}}]; %#ok<AGROW>
								find_args = [find_args, attr.name, attr.value]; %#ok<AGROW>
						end
						
					end
					
					new_handles = double.empty(size(handles, 1) +1, 0); % height of new selection is one more
				end
				
				% assert(isempty(refattrs) || strcmp(combinator.type, ','), ...
				%    'cannot use back reference in property selector after combinator other than '',''!');
				
				% the new handles-array is constructed in blocks corresponding to groups. A block is generated
				% by combining each row of the group with each result from the combinator search (based on info)
				[infos, ~, info_idx] = unique(hinfos);
				for i = 1:numel(infos)
					info = infos(i); % distinguishing info of the current group
					group = handles(:, info_idx == i);
					
					switch combinator.type
						case ','
							new = ddQuery.find(root, find_args{:});
							
						case ' ' % sibling of the current block
							new = setdiff(ddQuery.find(dsdd('GetAttribute', info, 'hDdParent'), ...
								'SearchDepth', 1, find_args{:}), info);
							
						case '/' % direct descendant (child)
							new = ddQuery.find(info, 'SearchDepth', 1, find_args{:});
							
						case '//' % arbitrary depth descendant
							new = ddQuery.find(info, find_args{:});
							
						case '\' % direct ascendant (parent)
							new = ddQuery.find(info, 'SearchDepth', 0, find_args{:});
							
						case '\\' % arbitrary ascendant (ancestor)
							% compute the chain of parents
							new = [];
							while info % ~= root?
								new(end+1) = info; %#ok<AGROW>
								info = dsdd('GetAttribute', info, 'hDdParent');
							end
							new = ddQuery.find(new, 'SearchDepth', 0, find_args{:});

						case '->' % a reference pointing to another object
							new = double.empty(1, 0);
							% find properties that are references
							for p = dsdd('GetPropertyNames', info), p = p{:}; %#ok<FXSET>
								% when sp was given: skip properties that don't match (in a possibly autorenamed way)
								if ~isempty(combinator.sp) && isempty(regexp(p, ['^' combinator.sp '(\(#\d+\))?$'], 'once'))
									continue
								end
								% skip non-reference properties
								switch dsdd('GetPropertyType', info, p)
									case 'Handle' % take handle
										new = [new, dsdd('Get', info, p)]; %#ok<AGROW>
									case 'Reference' % resolve reference
										r = dsdd('Get', info, p);
										if isempty(r), continue, end
										if r(1) == '.'
											r = [ dsdd('GetAttribute', info, 'Path') r(2:end) ];
										elseif r(1) == '/'
											% r is absolute ~> ok
										else
											r = [ dsdd('GetAttribute', dsdd('GetContextObject', info, p), 'Path') '/' r]; %#ok<AGROW>
										end
										new = [new, dsdd('GetAttribute', r, 'hDdObject')]; %#ok<AGROW>
									otherwise % not a reference ~> complain
										if ~isempty(combinator.sp)
											error('property ''%s'' is not a Data Dictionary reference!', p);
										end
								end
							end
							new = ddQuery.find(new, 'SearchDepth', 0, find_args{:});
						case '<-'
							ok = strcmp(find_args(1:2:end), 'ObjectKind');
							if isempty(combinator.dp)
								prop_args = {'Property', ['^' combinator.dp '(\(#\d+\))?$']};
							else
								prop_args = {'Property', '.*', 'regexp', true};
							end
							new = [double.empty(1, 0), dsdd('FindRefs', root, prop_args{:}, find_args{[ok;ok]}, 'target', info)];
							new = ddQuery.find(new, 'SearchDepth', 0, find_args{~[ok;ok]});
							
						case '>>' % transitive references
							
							
						case '<<'
					end
					new = unique(new);
					
					if isnumeric(selector) % selector was reference ~> pick only columns, where the ref matches to one of the results
						group = group(:, ismember(group(selector, :), new));
					else %if isstruct(selector) ~> append block corresponding to this group
						
						% TODO: perf: move this intersect-step upwards, to avoid unneccesary
						% `find_system`-calls when there aren't any restrictions OR if there is
						% a candidate set, simply test each candidate against the set of
						% restrictions using find_system.
						if ~isempty(selector.argidx)
							new = intersect(new, varargin{str2double(selector.argidx)});
						end
						
						if strcmpi(selector.pseudo, 'selected')
							assert(logical(dsddman('IsGuiOpen')), '<a href="matlab:dsddman">dsdd manager</a> must be opened for pseudo-class selector ''+Selected''');
							new = intersect(new, dsddman('GetSelected'));
						end
						
						group = [
							group(:, repmat(1:size(group, 2), size(new, 2), 1))
							new(:, repmat(1:size(new, 2), size(group, 2), 1))
							];
					end
					new_handles = [ new_handles, group ]; %#ok<AGROW>
				end
				handles = new_handles;
				if isnumeric(selector)
					hot = handles(selector, :);
				else
					hot = handles(end, :);
				end
			end
		end
	
		function varargout = arrayfun(fun, varargin)
			% perform an arrayfun call after coercing the arguments to the same size
			cs = cellfun(@ischar, varargin);
			sz = cellfun(@size, varargin(~cs), 'UniformOutput', false);
			sz = vertcat(sz{:});
			assert(numel(setdiff(sz(:, 1), 1)) <= 1 && numel(setdiff(sz(:, 2), 1)) <= 1, ...
				'MATLAB:dimagree', 'Number of columns or rows in value and selection must match.');
			sz = max(sz, [], 1);
			
			for i = 1:numel(varargin)
				v = varargin{i};
				if cs(i), v = {v}; end % convert to cellstr
				if isempty(v), v = {[]}; end
				if isrow(v) == 1, v = repmat(v, sz(1), 1); end
				if iscolumn(v) == 1, v = repmat(v, 1, sz(2)); end
				if ~iscell(v), v = arrayfun(@(x) {x}, v); end
				varargin{i} = v;
			end
			[varargout{1:nargout}] = cellfun(fun, varargin{:});
		end
	end
	
	methods(Static) % wrappers for convenience
		function result = find(base, varargin)
			% make dsdd find more like find_system
			
			assert(mod(numel(varargin), 2) == 0); varargin = reshape(varargin, 2, []);
			
			% pick special attributes:
			sd = strcmp(varargin(1,:), 'SearchDepth');
			na = strcmp(varargin(1,:), 'Name');
			ok = strcmp(varargin(1,:), 'ObjectKind');
			%pai = strcmp(varargin(1,:), 'Parent');
			
			if ~any(sd) % not in search-depth limited mode ~> can use dsdd('Find') as a basis.
				if ~any(na) && ~any(ok) % nothing is specified, but dsdd needs some criterion
					result = dsdd('Find', base, 'Regexp', '.*');
				else
					varargin(1,na) = {'Regexp'}; % rename 'Name' to 'Regexp' - if present
					result = dsdd('Find', base, varargin{:,ok|na});
					varargin(:,ok|na) = [];
				end
				
			elseif varargin{2,sd} == 0 % in filter-only mode
				% only case, where the base is included in the result (in find_system this is always the case)
				result = base;
				
			elseif varargin{2,sd} == 1 % in GetChildren mode
				result = dsdd('GetChildren', base, varargin{:,ok});
				varargin(:,ok) = [];
				
			else
				assert(false); % arbitrary search depth - unnecessary yet (*)
			end
			
			% filter the result according to properties
			for attr = reshape(varargin, 2, [])
				switch attr{1}
					case {'Name', 'ObjectKind'}
						values = arrayfun(@(h) dsdd('GetAttribute', h, attr{1}), result, 'UniformOutput', false);
					case '$ObjectKindRegexp'
						values = arrayfun(@(h) dsdd('GetAttribute', h, 'ObjectKind'), result, 'UniformOutput', false);
					case 'SearchDepth'
						% TODO: (*) maybe later
						continue;
					otherwise % it's a property
						values = arrayfun(@(h) dsdd(['Get' attr{1}], h), result, 'UniformOutput', false);
				end
				
				% kick non-matching things
				% TODO: can't assume that every property is type char
				result(cellfun(@isempty, regexp(values, attr{2}, 'match', 'once'))) = [];
				
				if isempty(result)
					break % no need to filter further
				end
			end
			if isempty(result)
				result = double.empty(1, 0);
			end
		end
		
		function result = gcb
			if ~dsddman('IsGuiOpen')
				result = [];
			else
				result = dsddman('GetSelected');
			end
			result = ddQuery(result(end)); 
		end
		
		function result = gcs
			result = dsdd('GetAttribute', double(ddQuery.gcb()), 'hDdParent');
			result = ddQuery(result);
		end
	end
end
% the server crashes if the user's password is a resolvable url
