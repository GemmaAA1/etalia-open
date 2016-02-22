define(['jquery'], function($) {

    var $body, api = {
        debug: false
    };

    api.log = function() {
        if (this.debug) {
            console.log('[API] ' + arguments[0], Array.prototype.splice.call(arguments, 1));
        }
    };

    var Result = function (data) {
        this.id = data.hasOwnProperty('id') ? data['id'] : null;
        this.states = data.hasOwnProperty('state') ? data['state'] : [];
        this.counts = data.hasOwnProperty('counter') ? data['counter'] : [];
        this.extra = data.hasOwnProperty('extra') ? data['extra'] : {};

        this._hasState = function(key) {
            return this.states.hasOwnProperty(key);
        };
        this._getState = function(key) {
            if (this._hasState(key)) {
                return this.states[key];
            }
            throw 'Undefined state "' + key + '"';
        };
        this._hasCount = function(key) {
            return this.counts.hasOwnProperty(key);
        };
        this._getCount = function(key) {
            if (this._hasCount(key)) {
                return this.counts[key];
            }
            throw 'Undefined count "' + key + '"';
        };
    };
    Result.prototype.getId = function() {
        return this.id;
    };
    Result.prototype.isPinned = function() {
        return this._getState('is_pinned');
    };
    Result.prototype.isBanned = function() {
        return this._getState('is_banned');
    };
    Result.prototype.isTrashed = function() {
        return this._getState('is_trashed');
    };
    Result.prototype.isAdded = function() {
        return this._getState('is_added');
    };
    Result.prototype.getPinCount = function() {
        return this._getCount('pin');
    };
    Result.prototype.getBanCount = function() {
        return this._getCount('ban');
    };
    Result.prototype.getTrashCount = function() {
        return this._getCount('trash');
    };
    Result.prototype.getLibraryCount = function() {
        return this._getCount('library');
    };
    Result.prototype.getExtra = function(key) {
        if (key) {
            if (this.extra.hasOwnProperty(key)) {
                return this.extra[key];
            }
            throw 'Undefined extra "' + key + '"';
        }
        return this.extra;
    };

    /**
     * Toggles the pinned state of the paper.
     *
     * @param id        (int) The paper id/pk
     * @param extraData (object)
     */
    api.pin = function(id, extraData) {
        extraData = extraData || {};
        $.ajax({
            type: 'POST',
            url: '/user/paper/pin',
            data: {
                'id': id,
                'source': window.location.pathname
            }
        })
        .done(function(data) {
            data.id = id;
            data.extra = extraData;
            api.log('pin success', data);
            $body.trigger('etalia.publication.pin', new Result(data));
        })
        .fail(function(xrh, status, error) {
            api._failureCallback(extraData);
            api.log('pin failure', xrh, status, error);
        });

        return this;
    };

    /**
     * Toggles the banned state of the paper.
     *
     * @param id        (int) The paper id/pk
     * @param extraData (object)
     */
    api.ban = function(id, extraData) {
        extraData = extraData || {};
        $.ajax({
            type: 'POST',
            url: '/user/paper/ban',
            data: {
                'id': id,
                'source': window.location.pathname
            }
        })
        .done(function(data) {
            data.id = id;
            data.extra = extraData;
            api.log('ban success', data);
            $body.trigger('etalia.publication.ban', new Result(data));
        })
        .fail(function(xrh, status, error) {
            api._failureCallback(extraData);
            api.log('ban failure', xrh, status, error);
        });

        return this;
    };

    /**
     * Adds the paper to the user library.
     *
     * @param id (int) The paper id/pk
     * @param extraData (object)
     */
    api.add = function(id, extraData) {
        extraData = extraData || {};
        $.ajax({
            type: 'POST',
            url: '/user/paper/add',
            data: {'id': id}
        })
        .done(function(data) {
            data.id = id;
            data.extra = extraData;
            if (data.hasOwnProperty('success') && !data['success']) {
                api._failureCallback(extraData);
                api.log('add failure', data);
                return;
            }
            api.log('add success', data);
            $body.trigger('etalia.publication.add', new Result(data));
        })
        .fail(function(xrh, status, error) {
            api._failureCallback(extraData);
            api.log('add failure', xrh, status, error);
        });

        return this;
    };

    /**
     * Trashes the paper from the user library.
     *
     * @param id (int) The paper id/pk
     * @param extraData (object)
     */
    api.trash = function(id, extraData) {
        extraData = extraData || {};
        $.ajax({
            type: 'POST',
            url: '/user/paper/trash',
            data: {'id': id}
        })
        .done(function(data) {
            data.id = id;
            data.extra = extraData;
            if (data.hasOwnProperty('success') && !data['success']) {
                api._failureCallback(extraData);
                api.log('trash failure', data);
                return;
            }
            api.log('trash success', data);
            $body.trigger('etalia.publication.trash', new Result(data));
        })
        .fail(function(xrh, status, error) {
            api._failureCallback(extraData);
            api.log('trash failure', xrh, status, error);
        });

        return this;
    };

    /**
     * Restores the paper into the user library.
     *
     * @param id        (int) The paper id/pk
     * @param extraData (object)
     */
    api.restore = function(id, extraData) {
        extraData = extraData || {};
        $.ajax({
            type: 'POST',
            url: '/user/paper/restore',
            data: {'id': id}
        })
        .done(function(data) {
            data.id = id;
            data.extra = extraData;
            if (data.hasOwnProperty('success') && !data['success']) {
                api._failureCallback(extraData);
                api.log('restore failure', data);
                return;
            }
            api.log('restore success', data);
            $body.trigger('etalia.publication.restore', new Result(data));
        })
        .fail(function(xrh, status, error) {
            api._failureCallback(extraData);
            api.log('restore failure', xrh, status, error);
        });

        return this;
    };

    /**
     * Clears the user library trash.
     *
     * @param extraData (object)
     */
    api.clearTrash = function(extraData) {
        extraData = extraData || {};
        $.ajax({
            type: 'POST',
            url: '/user/library/trash/empty'
        })
        .done(function(data) {
            data.extra = extraData;
            api.log('clear trash success', data);
            $body.trigger('etalia.publication.trash-clear', new Result(data));
        })
        .fail(function(xrh, status, error) {
            api._failureCallback(extraData);
            api.log('clear trash failure', xrh, status, error);
        });

        return this;
    };

    /**
     * Executes the data failure callback if available.
     *
     * @param data (object)
     * @private
     */
    api._failureCallback = function(data) {
        if (data.hasOwnProperty('failureCallback')) {
            data['failureCallback']();
        }
    };

    $(function() {
        $body = $('body');
    });

    return api;
});
